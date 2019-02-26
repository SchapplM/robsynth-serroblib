% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:37
% DurationCPUTime: 0.21s
% Computational Cost: add. (349->56), mult. (877->92), div. (0->0), fcn. (1158->12), ass. (0->41)
t30 = sin(pkin(11));
t35 = sin(qJ(2));
t39 = cos(qJ(2));
t50 = cos(pkin(11));
t43 = -t35 * t30 + t39 * t50;
t34 = sin(qJ(4));
t38 = cos(qJ(4));
t33 = sin(qJ(6));
t37 = cos(qJ(6));
t46 = t33 * r_i_i_C(1) + t37 * r_i_i_C(2) + qJ(5);
t49 = r_i_i_C(3) + pkin(10) + pkin(4);
t41 = t46 * t34 + t49 * t38 + pkin(3);
t56 = -pkin(5) - pkin(9);
t55 = t39 * pkin(2);
t32 = cos(pkin(6));
t54 = t32 * t39;
t31 = sin(pkin(6));
t36 = sin(qJ(1));
t52 = t36 * t31;
t40 = cos(qJ(1));
t51 = t40 * t31;
t24 = -t39 * t30 - t35 * t50;
t21 = t24 * t32;
t11 = -t40 * t21 + t36 * t43;
t47 = -t36 * t21 - t40 * t43;
t3 = t11 * t34 + t38 * t51;
t45 = -t11 * t38 + t34 * t51;
t44 = t37 * r_i_i_C(1) - t33 * r_i_i_C(2) - t56;
t42 = t43 * t32;
t29 = pkin(1) + t55;
t22 = t32 * t35 * pkin(2) + (-pkin(8) - qJ(3)) * t31;
t20 = t24 * t31;
t19 = t43 * t31;
t15 = -t20 * t34 - t32 * t38;
t13 = t40 * t24 - t36 * t42;
t10 = t36 * t24 + t40 * t42;
t8 = t34 * t52 - t38 * t47;
t7 = -t34 * t47 - t38 * t52;
t2 = -t13 * t37 + t7 * t33;
t1 = t13 * t33 + t7 * t37;
t4 = [-t11 * pkin(3) + t44 * t10 - t40 * t22 - t36 * t29 - t46 * t3 + t49 * t45 (-t40 * t35 - t36 * t54) * pkin(2) - t44 * t47 + t41 * t13, t52, t46 * t8 - t49 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -pkin(3) * t47 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t56 * t13 - t36 * t22 + t40 * t29 + t49 * t8 (-t36 * t35 + t40 * t54) * pkin(2) + t44 * t11 + t41 * t10, -t51, -t49 * t3 - t46 * t45, t3 (t10 * t33 + t3 * t37) * r_i_i_C(1) + (t10 * t37 - t3 * t33) * r_i_i_C(2); 0, t41 * t19 - t44 * t20 + t31 * t55, t32, t46 * (-t20 * t38 + t32 * t34) - t49 * t15, t15 (t15 * t37 + t19 * t33) * r_i_i_C(1) + (-t15 * t33 + t19 * t37) * r_i_i_C(2);];
Ja_transl  = t4;
