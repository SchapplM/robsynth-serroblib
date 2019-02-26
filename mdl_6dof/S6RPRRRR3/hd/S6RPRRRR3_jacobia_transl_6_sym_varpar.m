% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:15:59
% EndTime: 2019-02-26 21:16:00
% DurationCPUTime: 0.12s
% Computational Cost: add. (257->39), mult. (170->53), div. (0->0), fcn. (186->12), ass. (0->32)
t25 = cos(qJ(3));
t24 = sin(qJ(3));
t34 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t29 = t34 * t24;
t23 = qJ(4) + qJ(5);
t18 = cos(t23);
t11 = pkin(5) * t18 + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t39 = t25 * t9 + pkin(2) + t29;
t19 = qJ(6) + t23;
t13 = sin(t19);
t14 = cos(t19);
t21 = qJ(1) + pkin(11);
t16 = cos(t21);
t15 = sin(t21);
t33 = t15 * t25;
t5 = t13 * t33 + t16 * t14;
t6 = t16 * t13 - t14 * t33;
t38 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t32 = t16 * t25;
t7 = -t13 * t32 + t15 * t14;
t8 = t15 * t13 + t14 * t32;
t37 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t17 = sin(t23);
t36 = pkin(5) * t17;
t10 = t36 + sin(qJ(4)) * pkin(4);
t35 = pkin(7) + t10;
t31 = t17 * t25;
t28 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t27 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t26 = -t27 * t24 + t34 * t25;
t1 = [-sin(qJ(1)) * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t35 * t16 - t39 * t15, 0, t26 * t16, -t10 * t32 + t15 * t11 + t37 (t15 * t18 - t16 * t31) * pkin(5) + t37, t37; cos(qJ(1)) * pkin(1) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t35 * t15 + t39 * t16, 0, t26 * t15, -t10 * t33 - t16 * t11 + t38 (-t15 * t31 - t16 * t18) * pkin(5) + t38, t38; 0, 1, t27 * t25 + t29 (-t10 + t28) * t24 (t28 - t36) * t24, t28 * t24;];
Ja_transl  = t1;
