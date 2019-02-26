% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:53
% EndTime: 2019-02-26 22:35:53
% DurationCPUTime: 0.19s
% Computational Cost: add. (380->54), mult. (605->89), div. (0->0), fcn. (761->12), ass. (0->41)
t54 = pkin(4) + pkin(11) + r_i_i_C(3);
t36 = sin(qJ(6));
t40 = cos(qJ(6));
t50 = t36 * r_i_i_C(1) + t40 * r_i_i_C(2) + qJ(5);
t41 = cos(qJ(3));
t31 = t41 * pkin(3) + pkin(2);
t34 = qJ(3) + qJ(4);
t32 = sin(t34);
t33 = cos(t34);
t63 = t50 * t32 + t54 * t33 + t31;
t60 = pkin(5) + pkin(10) + pkin(9);
t35 = sin(pkin(6));
t38 = sin(qJ(2));
t59 = t35 * t38;
t39 = sin(qJ(1));
t58 = t35 * t39;
t42 = cos(qJ(2));
t57 = t35 * t42;
t43 = cos(qJ(1));
t56 = t35 * t43;
t55 = cos(pkin(6));
t52 = t43 * t55;
t24 = t38 * t52 + t39 * t42;
t12 = t24 * t33 - t32 * t56;
t53 = t39 * t55;
t37 = sin(qJ(3));
t51 = t35 * (pkin(3) * t37 + pkin(8));
t11 = t24 * t32 + t33 * t56;
t49 = t40 * r_i_i_C(1) - t36 * r_i_i_C(2) + t60;
t48 = -t54 * t11 + t50 * t12;
t26 = -t38 * t53 + t43 * t42;
t15 = t26 * t32 - t33 * t58;
t16 = t26 * t33 + t32 * t58;
t47 = -t54 * t15 + t50 * t16;
t21 = t32 * t59 - t55 * t33;
t46 = t50 * (t55 * t32 + t33 * t59) - t54 * t21;
t25 = t43 * t38 + t42 * t53;
t23 = t39 * t38 - t42 * t52;
t2 = t15 * t36 + t25 * t40;
t1 = t15 * t40 - t25 * t36;
t3 = [-t39 * pkin(1) - t50 * t11 - t54 * t12 - t49 * t23 - t24 * t31 + t43 * t51, -t25 * t63 + t49 * t26 (-t26 * t37 + t41 * t58) * pkin(3) + t47, t47, t15, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t43 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * qJ(5) + t54 * t16 + t60 * t25 + t26 * t31 + t39 * t51, -t23 * t63 + t49 * t24 (-t24 * t37 - t41 * t56) * pkin(3) + t48, t48, t11 (t11 * t40 - t23 * t36) * r_i_i_C(1) + (-t11 * t36 - t23 * t40) * r_i_i_C(2); 0 (t49 * t38 + t63 * t42) * t35 (-t37 * t59 + t55 * t41) * pkin(3) + t46, t46, t21 (t21 * t40 + t36 * t57) * r_i_i_C(1) + (-t21 * t36 + t40 * t57) * r_i_i_C(2);];
Ja_transl  = t3;
