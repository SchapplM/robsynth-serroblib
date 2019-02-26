% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:09
% EndTime: 2019-02-26 22:08:09
% DurationCPUTime: 0.27s
% Computational Cost: add. (350->69), mult. (908->117), div. (0->0), fcn. (1183->12), ass. (0->46)
t38 = sin(qJ(2));
t40 = cos(qJ(2));
t55 = cos(pkin(6));
t61 = cos(qJ(1));
t46 = t55 * t61;
t59 = sin(qJ(1));
t26 = t38 * t46 + t59 * t40;
t37 = sin(qJ(3));
t34 = sin(pkin(6));
t52 = t34 * t61;
t60 = cos(qJ(3));
t16 = t26 * t60 - t37 * t52;
t25 = t59 * t38 - t40 * t46;
t33 = sin(pkin(11));
t35 = cos(pkin(11));
t3 = t16 * t33 - t25 * t35;
t4 = t16 * t35 + t25 * t33;
t54 = -r_i_i_C(3) - pkin(10) + qJ(4);
t63 = t60 * pkin(3) + t54 * t37 + pkin(2);
t62 = pkin(4) + pkin(5);
t56 = t34 * t40;
t53 = t33 * t60;
t51 = t34 * t60;
t50 = t34 * t59;
t49 = t35 * t60;
t48 = t40 * t60;
t45 = t55 * t59;
t36 = sin(qJ(6));
t39 = cos(qJ(6));
t44 = t36 * r_i_i_C(1) + t39 * r_i_i_C(2) + qJ(5);
t43 = t39 * r_i_i_C(1) - t36 * r_i_i_C(2) + t62;
t15 = t26 * t37 + t61 * t51;
t41 = -t44 * t33 - t43 * t35 - pkin(3);
t28 = -t38 * t45 + t61 * t40;
t27 = t61 * t38 + t40 * t45;
t24 = t55 * t37 + t38 * t51;
t23 = t34 * t38 * t37 - t55 * t60;
t20 = t28 * t60 + t37 * t50;
t19 = t28 * t37 - t60 * t50;
t14 = t24 * t35 - t33 * t56;
t13 = t24 * t33 + t35 * t56;
t8 = t20 * t35 + t27 * t33;
t7 = t20 * t33 - t27 * t35;
t2 = t7 * t36 + t8 * t39;
t1 = -t8 * t36 + t7 * t39;
t5 = [-t59 * pkin(1) - t26 * pkin(2) - t16 * pkin(3) + pkin(8) * t52 - t25 * pkin(9) - t54 * t15 - t44 * t3 - t43 * t4, t28 * pkin(9) + t44 * (-t27 * t53 - t28 * t35) + t43 * (-t27 * t49 + t28 * t33) - t63 * t27, t41 * t19 + t54 * t20, t19, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t61 * pkin(1) + t28 * pkin(2) + t20 * pkin(3) + pkin(8) * t50 + t27 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t54 * t19 + t62 * t8, t26 * pkin(9) + t44 * (-t25 * t53 - t26 * t35) + t43 * (-t25 * t49 + t26 * t33) - t63 * t25, t41 * t15 + t54 * t16, t15, t3 (t3 * t39 - t4 * t36) * r_i_i_C(1) + (-t3 * t36 - t4 * t39) * r_i_i_C(2); 0 (t44 * (t33 * t48 - t35 * t38) + t43 * (t33 * t38 + t35 * t48) + t38 * pkin(9) + t63 * t40) * t34, t41 * t23 + t54 * t24, t23, t13 (t13 * t39 - t14 * t36) * r_i_i_C(1) + (-t13 * t36 - t14 * t39) * r_i_i_C(2);];
Ja_transl  = t5;
