% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRPR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:18
% EndTime: 2019-02-26 19:40:19
% DurationCPUTime: 0.18s
% Computational Cost: add. (355->48), mult. (889->87), div. (0->0), fcn. (1171->16), ass. (0->50)
t54 = sin(pkin(12));
t55 = sin(pkin(11));
t45 = t55 * t54;
t58 = cos(pkin(12));
t59 = cos(pkin(11));
t52 = t59 * t58;
t61 = cos(pkin(6));
t37 = -t61 * t52 + t45;
t56 = sin(pkin(7));
t57 = sin(pkin(6));
t49 = t57 * t56;
t60 = cos(pkin(7));
t66 = t37 * t60 + t59 * t49;
t47 = t55 * t58;
t50 = t59 * t54;
t38 = t61 * t47 + t50;
t46 = t55 * t57;
t65 = t38 * t60 - t56 * t46;
t64 = t58 * t60 * t57 + t61 * t56;
t63 = r_i_i_C(3) + pkin(10) + qJ(5);
t62 = cos(qJ(3));
t51 = t59 * t57;
t48 = t57 * t54;
t26 = pkin(13) + qJ(6);
t24 = sin(t26);
t25 = cos(t26);
t44 = -t25 * r_i_i_C(1) + t24 * r_i_i_C(2) - cos(pkin(13)) * pkin(5) - pkin(4);
t43 = sin(pkin(13)) * pkin(5) + t24 * r_i_i_C(1) + t25 * r_i_i_C(2) + pkin(9);
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t39 = -t63 * t29 + t44 * t31 - pkin(3);
t36 = -t58 * t49 + t61 * t60;
t33 = t38 * t56 + t60 * t46;
t32 = t37 * t56 - t60 * t51;
t30 = sin(qJ(3));
t19 = -t61 * t45 + t52;
t18 = t61 * t50 + t47;
t14 = t64 * t30 + t62 * t48;
t13 = t30 * t48 - t64 * t62;
t10 = t14 * t31 + t36 * t29;
t9 = t14 * t29 - t36 * t31;
t8 = t19 * t62 - t65 * t30;
t7 = t19 * t30 + t65 * t62;
t6 = t18 * t62 - t66 * t30;
t5 = t18 * t30 + t66 * t62;
t4 = t33 * t29 + t8 * t31;
t3 = t8 * t29 - t33 * t31;
t2 = t32 * t29 + t6 * t31;
t1 = t6 * t29 - t32 * t31;
t11 = [0, t46, t39 * t7 + t43 * t8, t44 * t3 + t63 * t4, t3 (-t4 * t24 + t7 * t25) * r_i_i_C(1) + (-t7 * t24 - t4 * t25) * r_i_i_C(2); 0, -t51, t39 * t5 + t43 * t6, t44 * t1 + t63 * t2, t1 (-t2 * t24 + t5 * t25) * r_i_i_C(1) + (-t2 * t25 - t5 * t24) * r_i_i_C(2); 1, t61, t39 * t13 + t43 * t14, t63 * t10 + t44 * t9, t9 (-t10 * t24 + t13 * t25) * r_i_i_C(1) + (-t10 * t25 - t13 * t24) * r_i_i_C(2);];
Ja_transl  = t11;
