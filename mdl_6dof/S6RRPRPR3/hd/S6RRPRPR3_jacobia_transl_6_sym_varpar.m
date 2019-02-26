% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:50
% EndTime: 2019-02-26 21:38:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (208->36), mult. (145->47), div. (0->0), fcn. (162->12), ass. (0->29)
t23 = qJ(2) + pkin(10);
t16 = sin(t23);
t37 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
t41 = cos(qJ(2)) * pkin(2) + t37 * t16;
t17 = cos(t23);
t22 = qJ(4) + pkin(11);
t11 = pkin(5) * cos(t22) + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t40 = t17 * t9 + pkin(1) + t41;
t18 = qJ(6) + t22;
t14 = cos(t18);
t27 = cos(qJ(1));
t13 = sin(t18);
t26 = sin(qJ(1));
t35 = t26 * t13;
t5 = t14 * t27 + t17 * t35;
t34 = t26 * t14;
t6 = t13 * t27 - t17 * t34;
t39 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t36 = t17 * t27;
t7 = -t13 * t36 + t34;
t8 = t14 * t36 + t35;
t38 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t10 = pkin(5) * sin(t22) + sin(qJ(4)) * pkin(4);
t33 = t10 + qJ(3) + pkin(7);
t30 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t29 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t28 = -sin(qJ(2)) * pkin(2) - t29 * t16 + t37 * t17;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t40 * t26 + t33 * t27, t28 * t27, t26, -t10 * t36 + t26 * t11 + t38, t27 * t16, t38; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t33 * t26 + t40 * t27, t28 * t26, -t27, -t26 * t17 * t10 - t11 * t27 + t39, t26 * t16, t39; 0, t17 * t29 + t41, 0 (-t10 + t30) * t16, -t17, t30 * t16;];
Ja_transl  = t1;
