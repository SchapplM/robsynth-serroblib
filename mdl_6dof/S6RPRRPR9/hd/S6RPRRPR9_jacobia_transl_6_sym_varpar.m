% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:28
% EndTime: 2019-02-26 21:05:29
% DurationCPUTime: 0.34s
% Computational Cost: add. (538->75), mult. (1252->124), div. (0->0), fcn. (1645->16), ass. (0->57)
t45 = cos(qJ(1));
t61 = sin(pkin(12));
t64 = cos(pkin(6));
t53 = t64 * t61;
t62 = cos(pkin(12));
t66 = sin(qJ(1));
t26 = t45 * t53 + t66 * t62;
t37 = sin(pkin(7));
t42 = sin(qJ(3));
t54 = t64 * t62;
t25 = -t45 * t54 + t66 * t61;
t63 = cos(pkin(7));
t58 = t25 * t63;
t38 = sin(pkin(6));
t67 = cos(qJ(3));
t60 = t38 * t67;
t11 = t45 * t37 * t60 + t26 * t42 + t67 * t58;
t65 = t45 * t38;
t12 = t26 * t67 + (-t37 * t65 - t58) * t42;
t57 = t38 * t63;
t21 = t25 * t37 - t45 * t57;
t36 = qJ(4) + pkin(13);
t34 = sin(t36);
t35 = cos(t36);
t4 = t12 * t35 + t21 * t34;
t40 = sin(qJ(6));
t43 = cos(qJ(6));
t76 = -t11 * t43 + t4 * t40;
t75 = -t11 * t40 - t4 * t43;
t41 = sin(qJ(4));
t74 = pkin(4) * t41 + pkin(9);
t71 = -t12 * t34 + t21 * t35;
t48 = t45 * t61 + t66 * t54;
t59 = t66 * t38;
t70 = -t37 * t59 + t48 * t63;
t69 = r_i_i_C(3) + pkin(11);
t56 = t64 * t37;
t52 = t63 * t62;
t51 = t43 * r_i_i_C(1) - t40 * r_i_i_C(2) + pkin(5);
t39 = -qJ(5) - pkin(10);
t50 = t40 * r_i_i_C(1) + t43 * r_i_i_C(2) - t39;
t44 = cos(qJ(4));
t33 = t44 * pkin(4) + pkin(3);
t49 = -t69 * t34 - t51 * t35 - t33;
t46 = t48 * t37 + t66 * t57;
t27 = t45 * t62 - t66 * t53;
t24 = -t38 * t62 * t37 + t64 * t63;
t19 = t42 * t56 + (t42 * t52 + t67 * t61) * t38;
t18 = t38 * t61 * t42 - t52 * t60 - t67 * t56;
t16 = t27 * t67 - t70 * t42;
t15 = t27 * t42 + t70 * t67;
t10 = t19 * t35 + t24 * t34;
t8 = t16 * t35 + t34 * t46;
t7 = t16 * t34 - t35 * t46;
t2 = t15 * t40 + t8 * t43;
t1 = t15 * t43 - t8 * t40;
t3 = [-t66 * pkin(1) - t26 * pkin(2) - t4 * pkin(5) + t75 * r_i_i_C(1) + t76 * r_i_i_C(2) + qJ(2) * t65 + t11 * t39 - t12 * t33 - t74 * t21 + t69 * t71, t59, t49 * t15 + t50 * t16, t69 * t8 + (-t16 * t41 + t44 * t46) * pkin(4) - t51 * t7, t15, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t45 * pkin(1) + t27 * pkin(2) + t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t59 - t15 * t39 + t16 * t33 + t74 * t46 + t69 * t7, -t65, t49 * t11 + t50 * t12, t69 * t4 + (-t12 * t41 + t21 * t44) * pkin(4) + t51 * t71, t11, -t76 * r_i_i_C(1) + t75 * r_i_i_C(2); 0, t64, t49 * t18 + t50 * t19, t69 * t10 + (-t19 * t41 + t24 * t44) * pkin(4) + t51 * (-t19 * t34 + t24 * t35) t18 (-t10 * t40 + t18 * t43) * r_i_i_C(1) + (-t10 * t43 - t18 * t40) * r_i_i_C(2);];
Ja_transl  = t3;
