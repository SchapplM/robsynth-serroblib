% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:45
% EndTime: 2019-02-26 21:55:46
% DurationCPUTime: 0.17s
% Computational Cost: add. (214->50), mult. (420->80), div. (0->0), fcn. (539->12), ass. (0->39)
t28 = qJ(4) + qJ(5);
t26 = sin(t28);
t30 = sin(pkin(6));
t38 = cos(qJ(1));
t44 = t38 * t30;
t22 = t26 * t44;
t27 = cos(t28);
t32 = cos(pkin(6));
t29 = sin(pkin(12));
t31 = cos(pkin(12));
t34 = sin(qJ(2));
t37 = cos(qJ(2));
t41 = t37 * t29 + t34 * t31;
t18 = t41 * t32;
t20 = t34 * t29 - t37 * t31;
t35 = sin(qJ(1));
t9 = t38 * t18 - t35 * t20;
t51 = (-t9 * t26 - t27 * t44) * r_i_i_C(1) + (-t9 * t27 + t22) * r_i_i_C(2);
t42 = t35 * t18 + t38 * t20;
t45 = t35 * t30;
t5 = t26 * t42 + t27 * t45;
t6 = t26 * t45 - t27 * t42;
t50 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t33 = sin(qJ(4));
t49 = pkin(4) * t33;
t48 = t37 * pkin(2);
t47 = r_i_i_C(3) + pkin(10) + pkin(9);
t46 = t32 * t37;
t16 = t41 * t30;
t43 = (-t16 * t26 + t32 * t27) * r_i_i_C(1) + (-t16 * t27 - t32 * t26) * r_i_i_C(2);
t36 = cos(qJ(4));
t24 = t36 * pkin(4) + pkin(3);
t40 = t27 * r_i_i_C(1) - t26 * r_i_i_C(2) + t24;
t25 = pkin(1) + t48;
t19 = t32 * t34 * pkin(2) + (-pkin(8) - qJ(3)) * t30;
t17 = t20 * t32;
t11 = t35 * t17 - t38 * t41;
t8 = -t38 * t17 - t35 * t41;
t1 = [t22 * r_i_i_C(1) - t35 * t25 - t40 * t9 + t47 * t8 + (-t19 + (r_i_i_C(2) * t27 + t49) * t30) * t38, -t47 * t42 + (-t34 * t38 - t35 * t46) * pkin(2) + t40 * t11, t45 (t33 * t42 + t36 * t45) * pkin(4) + t50, t50, 0; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t42 * t24 + t38 * t25 + (t30 * t49 - t19) * t35 - t47 * t11, t47 * t9 + (-t34 * t35 + t38 * t46) * pkin(2) + t40 * t8, -t44 (-t33 * t9 - t36 * t44) * pkin(4) + t51, t51, 0; 0, t47 * t16 + (-t20 * t40 + t48) * t30, t32 (-t16 * t33 + t32 * t36) * pkin(4) + t43, t43, 0;];
Ja_transl  = t1;
