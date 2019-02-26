% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:30
% EndTime: 2019-02-26 20:54:31
% DurationCPUTime: 0.32s
% Computational Cost: add. (517->68), mult. (1194->112), div. (0->0), fcn. (1572->16), ass. (0->55)
t43 = cos(qJ(1));
t60 = sin(pkin(12));
t64 = cos(pkin(6));
t52 = t64 * t60;
t62 = cos(pkin(12));
t66 = sin(qJ(1));
t26 = t43 * t52 + t66 * t62;
t41 = sin(qJ(3));
t67 = cos(qJ(3));
t54 = t64 * t62;
t25 = -t43 * t54 + t66 * t60;
t38 = sin(pkin(6));
t61 = sin(pkin(7));
t56 = t38 * t61;
t63 = cos(pkin(7));
t72 = t25 * t63 + t43 * t56;
t11 = t26 * t41 + t72 * t67;
t12 = t26 * t67 - t72 * t41;
t57 = t38 * t63;
t21 = t25 * t61 - t43 * t57;
t36 = pkin(13) + qJ(5);
t34 = sin(t36);
t35 = cos(t36);
t4 = t12 * t35 + t21 * t34;
t40 = sin(qJ(6));
t42 = cos(qJ(6));
t77 = -t11 * t42 + t4 * t40;
t76 = -t11 * t40 - t4 * t42;
t75 = pkin(9) + pkin(4) * sin(pkin(13));
t71 = -t12 * t34 + t21 * t35;
t46 = t43 * t60 + t66 * t54;
t70 = t46 * t63 - t66 * t56;
t69 = r_i_i_C(3) + pkin(11);
t65 = t43 * t38;
t59 = t66 * t38;
t53 = t64 * t61;
t51 = t63 * t62;
t50 = t42 * r_i_i_C(1) - t40 * r_i_i_C(2) + pkin(5);
t39 = -pkin(10) - qJ(4);
t49 = t40 * r_i_i_C(1) + t42 * r_i_i_C(2) - t39;
t33 = cos(pkin(13)) * pkin(4) + pkin(3);
t47 = -t69 * t34 - t50 * t35 - t33;
t44 = t46 * t61 + t66 * t57;
t27 = t43 * t62 - t66 * t52;
t24 = -t62 * t56 + t64 * t63;
t19 = t41 * t53 + (t41 * t51 + t67 * t60) * t38;
t18 = -t67 * t53 + (t41 * t60 - t51 * t67) * t38;
t16 = t27 * t67 - t70 * t41;
t15 = t27 * t41 + t70 * t67;
t10 = t19 * t35 + t24 * t34;
t8 = t16 * t35 + t34 * t44;
t7 = t16 * t34 - t35 * t44;
t2 = t15 * t40 + t8 * t42;
t1 = t15 * t42 - t8 * t40;
t3 = [-t66 * pkin(1) - t26 * pkin(2) - t4 * pkin(5) + t76 * r_i_i_C(1) + t77 * r_i_i_C(2) + qJ(2) * t65 + t11 * t39 - t12 * t33 - t75 * t21 + t69 * t71, t59, t47 * t15 + t49 * t16, t15, -t50 * t7 + t69 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t43 * pkin(1) + t27 * pkin(2) + t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t59 - t15 * t39 + t16 * t33 + t75 * t44 + t69 * t7, -t65, t47 * t11 + t49 * t12, t11, t69 * t4 + t50 * t71, -t77 * r_i_i_C(1) + t76 * r_i_i_C(2); 0, t64, t47 * t18 + t49 * t19, t18, t69 * t10 + t50 * (-t19 * t34 + t24 * t35) (-t10 * t40 + t18 * t42) * r_i_i_C(1) + (-t10 * t42 - t18 * t40) * r_i_i_C(2);];
Ja_transl  = t3;
