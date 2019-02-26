% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:27
% EndTime: 2019-02-26 20:06:27
% DurationCPUTime: 0.18s
% Computational Cost: add. (274->44), mult. (408->76), div. (0->0), fcn. (514->14), ass. (0->38)
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t28 = pkin(12) + qJ(5);
t26 = qJ(6) + t28;
t22 = sin(t26);
t23 = cos(t26);
t25 = cos(t28);
t37 = r_i_i_C(1) * t23 - r_i_i_C(2) * t22 + pkin(5) * t25 + cos(pkin(12)) * pkin(4) + pkin(3);
t46 = r_i_i_C(3) + pkin(10) + pkin(9) + qJ(4);
t50 = t46 * t31 + t37 * t33 + pkin(2);
t29 = sin(pkin(11));
t32 = sin(qJ(2));
t34 = cos(qJ(2));
t41 = cos(pkin(11));
t42 = cos(pkin(6));
t38 = t42 * t41;
t11 = t29 * t32 - t34 * t38;
t12 = t29 * t34 + t32 * t38;
t30 = sin(pkin(6));
t39 = t30 * t41;
t8 = t12 * t33 - t31 * t39;
t49 = (t11 * t23 - t8 * t22) * r_i_i_C(1) + (-t11 * t22 - t8 * t23) * r_i_i_C(2);
t40 = t29 * t42;
t14 = -t32 * t40 + t41 * t34;
t45 = t30 * t31;
t10 = t14 * t33 + t29 * t45;
t13 = t41 * t32 + t34 * t40;
t48 = (-t10 * t22 + t13 * t23) * r_i_i_C(1) + (-t10 * t23 - t13 * t22) * r_i_i_C(2);
t44 = t30 * t33;
t16 = t42 * t31 + t32 * t44;
t43 = t30 * t34;
t47 = (-t16 * t22 - t23 * t43) * r_i_i_C(1) + (-t16 * t23 + t22 * t43) * r_i_i_C(2);
t24 = sin(t28);
t36 = t22 * r_i_i_C(1) + t23 * r_i_i_C(2) + pkin(8) + pkin(5) * t24 + sin(pkin(12)) * pkin(4);
t15 = t32 * t45 - t42 * t33;
t9 = t14 * t31 - t29 * t44;
t7 = t12 * t31 + t33 * t39;
t1 = [0, -t13 * t50 + t36 * t14, t46 * t10 - t37 * t9, t9 (-t10 * t24 + t13 * t25) * pkin(5) + t48, t48; 0, -t11 * t50 + t36 * t12, -t37 * t7 + t46 * t8, t7 (t11 * t25 - t24 * t8) * pkin(5) + t49, t49; 1 (t36 * t32 + t34 * t50) * t30, -t37 * t15 + t46 * t16, t15 (-t16 * t24 - t25 * t43) * pkin(5) + t47, t47;];
Ja_transl  = t1;
