% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:10
% EndTime: 2019-02-26 21:55:10
% DurationCPUTime: 0.16s
% Computational Cost: add. (250->40), mult. (175->52), div. (0->0), fcn. (193->12), ass. (0->35)
t23 = qJ(2) + pkin(11);
t16 = sin(t23);
t42 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t47 = cos(qJ(2)) * pkin(2) + t42 * t16;
t17 = cos(t23);
t25 = qJ(4) + qJ(5);
t19 = cos(t25);
t11 = pkin(5) * t19 + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t46 = t17 * t9 + pkin(1) + t47;
t20 = qJ(6) + t25;
t14 = cos(t20);
t29 = cos(qJ(1));
t36 = t29 * t14;
t13 = sin(t20);
t28 = sin(qJ(1));
t39 = t28 * t13;
t5 = t17 * t39 + t36;
t37 = t29 * t13;
t38 = t28 * t14;
t6 = -t17 * t38 + t37;
t45 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t17 * t37 + t38;
t8 = t17 * t36 + t39;
t44 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t18 = sin(t25);
t43 = pkin(5) * t18;
t41 = t17 * t28;
t40 = t17 * t29;
t10 = t43 + sin(qJ(4)) * pkin(4);
t35 = t10 + qJ(3) + pkin(7);
t32 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t31 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t30 = -sin(qJ(2)) * pkin(2) - t31 * t16 + t42 * t17;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t46 * t28 + t35 * t29, t30 * t29, t28, -t10 * t40 + t28 * t11 + t44 (-t18 * t40 + t19 * t28) * pkin(5) + t44, t44; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t35 * t28 + t29 * t46, t30 * t28, -t29, -t10 * t41 - t29 * t11 + t45 (-t18 * t41 - t19 * t29) * pkin(5) + t45, t45; 0, t31 * t17 + t47, 0 (-t10 + t32) * t16 (t32 - t43) * t16, t32 * t16;];
Ja_transl  = t1;
