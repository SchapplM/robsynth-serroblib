% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:42:33
% EndTime: 2019-02-26 21:42:33
% DurationCPUTime: 0.15s
% Computational Cost: add. (218->46), mult. (352->72), div. (0->0), fcn. (446->12), ass. (0->33)
t17 = cos(pkin(11)) * pkin(3) + pkin(2);
t20 = pkin(11) + qJ(4);
t18 = sin(t20);
t19 = cos(t20);
t21 = sin(pkin(12));
t24 = cos(pkin(12));
t32 = r_i_i_C(1) * t24 - r_i_i_C(2) * t21 + pkin(4);
t38 = r_i_i_C(3) + qJ(5);
t42 = t38 * t18 + t32 * t19 + t17;
t23 = sin(pkin(6));
t26 = sin(qJ(2));
t41 = t23 * t26;
t27 = sin(qJ(1));
t40 = t23 * t27;
t29 = cos(qJ(1));
t39 = t23 * t29;
t37 = cos(pkin(6));
t28 = cos(qJ(2));
t34 = t29 * t37;
t10 = t26 * t34 + t27 * t28;
t36 = t10 * t19 - t18 * t39;
t35 = t27 * t37;
t33 = t23 * (pkin(3) * sin(pkin(11)) + pkin(8));
t25 = -pkin(9) - qJ(3);
t31 = t21 * r_i_i_C(1) + t24 * r_i_i_C(2) - t25;
t1 = t10 * t18 + t19 * t39;
t12 = -t26 * t35 + t29 * t28;
t11 = t29 * t26 + t28 * t35;
t9 = t27 * t26 - t28 * t34;
t7 = t18 * t41 - t37 * t19;
t6 = t12 * t19 + t18 * t40;
t5 = t12 * t18 - t19 * t40;
t2 = [(-t9 * t21 - t24 * t36) * r_i_i_C(1) + (t21 * t36 - t9 * t24) * r_i_i_C(2) - t36 * pkin(4) - t10 * t17 + t9 * t25 - t27 * pkin(1) - t38 * t1 + t29 * t33, -t11 * t42 + t31 * t12, t11, -t32 * t5 + t38 * t6, t5, 0; (t11 * t21 + t6 * t24) * r_i_i_C(1) + (t11 * t24 - t6 * t21) * r_i_i_C(2) + t6 * pkin(4) + t12 * t17 - t11 * t25 + t29 * pkin(1) + t38 * t5 + t27 * t33, t31 * t10 - t42 * t9, t9, -t32 * t1 + t38 * t36, t1, 0; 0 (t31 * t26 + t42 * t28) * t23, -t23 * t28, t38 * (t37 * t18 + t19 * t41) - t32 * t7, t7, 0;];
Ja_transl  = t2;
