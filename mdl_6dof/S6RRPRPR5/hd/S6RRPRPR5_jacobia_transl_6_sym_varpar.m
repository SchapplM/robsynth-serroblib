% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR5
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
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:07
% EndTime: 2019-02-26 21:40:08
% DurationCPUTime: 0.22s
% Computational Cost: add. (370->58), mult. (809->94), div. (0->0), fcn. (1065->14), ass. (0->44)
t36 = sin(pkin(11));
t41 = sin(qJ(2));
t44 = cos(qJ(2));
t54 = cos(pkin(11));
t49 = -t41 * t36 + t44 * t54;
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t30 = cos(pkin(12)) * pkin(5) + pkin(4);
t34 = pkin(12) + qJ(6);
t32 = sin(t34);
t33 = cos(t34);
t50 = t33 * r_i_i_C(1) - t32 * r_i_i_C(2) + t30;
t59 = r_i_i_C(3) + pkin(10) + qJ(5);
t46 = t59 * t40 + t50 * t43 + pkin(3);
t60 = t44 * pkin(2);
t38 = cos(pkin(6));
t58 = t38 * t44;
t37 = sin(pkin(6));
t42 = sin(qJ(1));
t56 = t42 * t37;
t45 = cos(qJ(1));
t55 = t45 * t37;
t53 = -sin(pkin(12)) * pkin(5) - pkin(9);
t24 = -t44 * t36 - t41 * t54;
t21 = t24 * t38;
t11 = -t45 * t21 + t42 * t49;
t4 = t11 * t43 - t40 * t55;
t51 = -t42 * t21 - t45 * t49;
t3 = t11 * t40 + t43 * t55;
t48 = t32 * r_i_i_C(1) + t33 * r_i_i_C(2) - t53;
t47 = t49 * t38;
t31 = pkin(1) + t60;
t22 = t38 * t41 * pkin(2) + (-pkin(8) - qJ(3)) * t37;
t20 = t24 * t37;
t19 = t49 * t37;
t16 = -t20 * t43 + t38 * t40;
t15 = -t20 * t40 - t38 * t43;
t13 = t45 * t24 - t42 * t47;
t10 = t42 * t24 + t45 * t47;
t8 = t40 * t56 - t43 * t51;
t7 = -t40 * t51 - t43 * t56;
t2 = -t13 * t32 + t8 * t33;
t1 = -t13 * t33 - t8 * t32;
t5 = [-t11 * pkin(3) + t48 * t10 - t45 * t22 - t59 * t3 - t42 * t31 - t50 * t4 (-t45 * t41 - t42 * t58) * pkin(2) - t48 * t51 + t46 * t13, t56, -t50 * t7 + t59 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -pkin(3) * t51 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t53 * t13 - t42 * t22 + t8 * t30 + t45 * t31 + t59 * t7 (-t42 * t41 + t45 * t58) * pkin(2) + t48 * t11 + t46 * t10, -t55, -t50 * t3 + t59 * t4, t3 (-t10 * t33 - t4 * t32) * r_i_i_C(1) + (t10 * t32 - t4 * t33) * r_i_i_C(2); 0, t46 * t19 - t48 * t20 + t37 * t60, t38, -t50 * t15 + t59 * t16, t15 (-t16 * t32 - t19 * t33) * r_i_i_C(1) + (-t16 * t33 + t19 * t32) * r_i_i_C(2);];
Ja_transl  = t5;
