% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR13_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR13_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:03
% EndTime: 2019-02-26 22:23:04
% DurationCPUTime: 0.22s
% Computational Cost: add. (317->68), mult. (734->116), div. (0->0), fcn. (944->14), ass. (0->49)
t68 = pkin(10) + pkin(4) * sin(pkin(13));
t33 = cos(pkin(13)) * pkin(4) + pkin(3);
t36 = pkin(13) + qJ(5);
t34 = sin(t36);
t35 = cos(t36);
t49 = t35 * r_i_i_C(1) - t34 * r_i_i_C(2) + t33;
t67 = t34 * r_i_i_C(1) + t35 * r_i_i_C(2) + t68;
t65 = r_i_i_C(3) + pkin(11) + qJ(4);
t43 = sin(qJ(2));
t44 = sin(qJ(1));
t46 = cos(qJ(2));
t47 = cos(qJ(1));
t55 = cos(pkin(6));
t51 = t47 * t55;
t24 = t43 * t51 + t44 * t46;
t42 = sin(qJ(3));
t64 = t24 * t42;
t39 = sin(pkin(6));
t63 = t39 * t44;
t62 = t39 * t47;
t40 = cos(pkin(7));
t61 = t40 * t42;
t45 = cos(qJ(3));
t60 = t40 * t45;
t59 = t42 * t43;
t58 = t42 * t46;
t57 = t43 * t45;
t56 = t45 * t46;
t38 = sin(pkin(7));
t54 = t38 * t63;
t53 = t38 * t62;
t52 = t44 * t55;
t50 = t55 * t38;
t23 = t44 * t43 - t46 * t51;
t15 = -t23 * t38 + t40 * t62;
t25 = -t47 * t43 - t46 * t52;
t17 = -t25 * t38 + t40 * t63;
t6 = t23 * t61 - t24 * t45 + t42 * t53;
t48 = t67 * t38;
t26 = -t43 * t52 + t47 * t46;
t22 = -t39 * t46 * t38 + t55 * t40;
t14 = t42 * t50 + (t40 * t58 + t57) * t39;
t13 = -t45 * t50 + (-t40 * t56 + t59) * t39;
t8 = t26 * t45 + (t25 * t40 + t54) * t42;
t7 = -t25 * t60 + t26 * t42 - t45 * t54;
t3 = t23 * t60 + t45 * t53 + t64;
t2 = t17 * t34 + t8 * t35;
t1 = t17 * t35 - t8 * t34;
t4 = [-t24 * pkin(2) - t44 * pkin(1) + pkin(9) * t62 + t65 * (-t64 + (-t23 * t40 - t53) * t45) + t49 * t6 + t67 * t15, t25 * pkin(2) + t49 * (t25 * t45 - t26 * t61) + t26 * t48 + t65 * (t25 * t42 + t26 * t60) -t49 * t7 + t65 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t47 * pkin(1) + t26 * pkin(2) + pkin(9) * t63 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t68 * t17 + t8 * t33 + t65 * t7, -t23 * pkin(2) + t65 * (-t23 * t42 + t24 * t60) + t49 * (-t23 * t45 - t24 * t61) + t24 * t48, -t49 * t3 - t6 * t65, t3 (-t15 * t35 + t6 * t34) * r_i_i_C(1) + (t15 * t34 + t6 * t35) * r_i_i_C(2), 0; 0 (t65 * (t40 * t57 + t58) + t49 * (-t40 * t59 + t56) + t46 * pkin(2) + t43 * t48) * t39, -t49 * t13 + t65 * t14, t13 (-t14 * t34 + t22 * t35) * r_i_i_C(1) + (-t14 * t35 - t22 * t34) * r_i_i_C(2), 0;];
Ja_transl  = t4;
