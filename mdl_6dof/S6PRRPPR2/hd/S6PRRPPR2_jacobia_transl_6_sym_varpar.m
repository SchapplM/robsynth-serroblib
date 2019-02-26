% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:45
% EndTime: 2019-02-26 19:58:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (214->43), mult. (356->73), div. (0->0), fcn. (449->12), ass. (0->31)
t17 = qJ(3) + pkin(11);
t15 = sin(t17);
t16 = cos(t17);
t25 = cos(qJ(3));
t21 = sin(qJ(6));
t24 = cos(qJ(6));
t29 = t21 * r_i_i_C(1) + t24 * r_i_i_C(2) + qJ(5);
t33 = pkin(4) + pkin(9) + r_i_i_C(3);
t39 = t25 * pkin(3) + t29 * t15 + t33 * t16 + pkin(2);
t18 = sin(pkin(10));
t19 = sin(pkin(6));
t38 = t18 * t19;
t23 = sin(qJ(2));
t37 = t19 * t23;
t26 = cos(qJ(2));
t36 = t19 * t26;
t35 = cos(pkin(6));
t34 = cos(pkin(10));
t32 = t18 * t35;
t31 = t19 * t34;
t30 = t35 * t34;
t28 = t24 * r_i_i_C(1) - t21 * r_i_i_C(2) + pkin(5) + pkin(8) + qJ(4);
t22 = sin(qJ(3));
t10 = -t23 * t32 + t34 * t26;
t9 = t34 * t23 + t26 * t32;
t8 = t18 * t26 + t23 * t30;
t7 = t18 * t23 - t26 * t30;
t5 = t15 * t37 - t35 * t16;
t3 = t10 * t15 - t16 * t38;
t1 = t8 * t15 + t16 * t31;
t2 = [0, t28 * t10 - t39 * t9 (-t10 * t22 + t25 * t38) * pkin(3) + t29 * (t10 * t16 + t15 * t38) - t33 * t3, t9, t3 (-t9 * t21 + t3 * t24) * r_i_i_C(1) + (-t3 * t21 - t9 * t24) * r_i_i_C(2); 0, t28 * t8 - t39 * t7 (-t22 * t8 - t25 * t31) * pkin(3) + t29 * (-t15 * t31 + t8 * t16) - t33 * t1, t7, t1 (t1 * t24 - t7 * t21) * r_i_i_C(1) + (-t1 * t21 - t7 * t24) * r_i_i_C(2); 1 (t28 * t23 + t39 * t26) * t19 (-t22 * t37 + t35 * t25) * pkin(3) + t29 * (t35 * t15 + t16 * t37) - t33 * t5, -t36, t5 (t21 * t36 + t5 * t24) * r_i_i_C(1) + (-t5 * t21 + t24 * t36) * r_i_i_C(2);];
Ja_transl  = t2;
