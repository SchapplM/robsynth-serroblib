% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:03
% EndTime: 2019-02-26 19:47:03
% DurationCPUTime: 0.18s
% Computational Cost: add. (260->43), mult. (577->77), div. (0->0), fcn. (759->14), ass. (0->37)
t29 = sin(pkin(11));
t36 = sin(qJ(2));
t38 = cos(qJ(2));
t46 = cos(pkin(11));
t42 = -t36 * t29 + t38 * t46;
t35 = sin(qJ(4));
t37 = cos(qJ(4));
t27 = pkin(12) + qJ(6);
t25 = sin(t27);
t26 = cos(t27);
t43 = t26 * r_i_i_C(1) - t25 * r_i_i_C(2) + cos(pkin(12)) * pkin(5) + pkin(4);
t51 = r_i_i_C(3) + pkin(9) + qJ(5);
t39 = t51 * t35 + t43 * t37 + pkin(3);
t30 = sin(pkin(10));
t31 = sin(pkin(6));
t50 = t30 * t31;
t32 = cos(pkin(10));
t49 = t32 * t31;
t33 = cos(pkin(6));
t48 = t33 * t38;
t19 = -t38 * t29 - t36 * t46;
t17 = t19 * t33;
t7 = -t32 * t17 + t30 * t42;
t44 = -t30 * t17 - t32 * t42;
t41 = sin(pkin(12)) * pkin(5) + t25 * r_i_i_C(1) + t26 * r_i_i_C(2) + pkin(8);
t40 = t42 * t33;
t16 = t19 * t31;
t15 = t42 * t31;
t12 = -t16 * t37 + t33 * t35;
t11 = -t16 * t35 - t33 * t37;
t9 = t32 * t19 - t30 * t40;
t6 = t30 * t19 + t32 * t40;
t4 = t35 * t50 - t37 * t44;
t3 = -t35 * t44 - t37 * t50;
t2 = -t35 * t49 + t7 * t37;
t1 = t7 * t35 + t37 * t49;
t5 = [0 (-t30 * t48 - t32 * t36) * pkin(2) - t41 * t44 + t39 * t9, t50, -t43 * t3 + t51 * t4, t3 (-t4 * t25 - t9 * t26) * r_i_i_C(1) + (t9 * t25 - t4 * t26) * r_i_i_C(2); 0 (-t30 * t36 + t32 * t48) * pkin(2) + t41 * t7 + t39 * t6, -t49, -t43 * t1 + t51 * t2, t1 (-t2 * t25 - t6 * t26) * r_i_i_C(1) + (-t2 * t26 + t6 * t25) * r_i_i_C(2); 1, t31 * t38 * pkin(2) + t39 * t15 - t41 * t16, t33, -t43 * t11 + t51 * t12, t11 (-t12 * t25 - t15 * t26) * r_i_i_C(1) + (-t12 * t26 + t15 * t25) * r_i_i_C(2);];
Ja_transl  = t5;
