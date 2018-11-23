% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10_jacobia_transl_3_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_transl_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobia_transl_3_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_transl_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:16
% EndTime: 2018-11-23 11:27:17
% DurationCPUTime: 0.17s
% Computational Cost: add. (312->59), mult. (325->85), div. (0->0), fcn. (316->18), ass. (0->43)
t54 = pkin(11) + r_i_i_C(3);
t24 = sin(pkin(6));
t29 = sin(qJ(1));
t53 = t24 * t29;
t32 = cos(qJ(1));
t52 = t24 * t32;
t51 = pkin(6) - qJ(2);
t50 = pkin(6) + qJ(2);
t49 = pkin(7) - qJ(3);
t48 = pkin(7) + qJ(3);
t23 = sin(pkin(7));
t47 = t54 * t23;
t46 = cos(t51);
t45 = cos(t49);
t44 = sin(t51);
t43 = sin(t49);
t42 = cos(t50) / 0.2e1;
t41 = cos(t48) / 0.2e1;
t40 = sin(t50) / 0.2e1;
t39 = sin(t48) / 0.2e1;
t12 = t40 - t44 / 0.2e1;
t31 = cos(qJ(2));
t3 = t32 * t12 + t29 * t31;
t7 = t29 * t12 - t32 * t31;
t14 = t45 / 0.2e1 + t41;
t28 = sin(qJ(2));
t33 = t46 / 0.2e1 + t42;
t2 = t29 * t28 - t32 * t33;
t27 = sin(qJ(3));
t9 = t39 + t43 / 0.2e1;
t36 = t2 * t14 + t3 * t27 + t9 * t52;
t10 = t39 - t43 / 0.2e1;
t13 = t41 - t45 / 0.2e1;
t30 = cos(qJ(3));
t35 = t2 * t10 - t13 * t52 - t3 * t30;
t5 = -t32 * t28 - t29 * t33;
t34 = -t5 * t10 + t13 * t53 + t30 * t7;
t26 = cos(pkin(6));
t25 = cos(pkin(7));
t15 = t42 - t46 / 0.2e1;
t11 = t40 + t44 / 0.2e1;
t1 = t5 * t14 + t27 * t7 + t9 * t53;
t4 = [-t29 * pkin(1) - t3 * pkin(2) + pkin(10) * t52 + t35 * r_i_i_C(1) + t36 * r_i_i_C(2) + t54 * (-t2 * t23 + t25 * t52) (t7 * t10 + t5 * t30) * r_i_i_C(1) + (t7 * t14 - t5 * t27) * r_i_i_C(2) + t5 * pkin(2) - t7 * t47, t1 * r_i_i_C(1) + t34 * r_i_i_C(2), 0, 0, 0; t32 * pkin(1) - t7 * pkin(2) + pkin(10) * t53 - t34 * r_i_i_C(1) + t1 * r_i_i_C(2) + t54 * (-t5 * t23 + t25 * t53) (-t10 * t3 - t2 * t30) * r_i_i_C(1) + (-t14 * t3 + t2 * t27) * r_i_i_C(2) - t2 * pkin(2) + t3 * t47, -t36 * r_i_i_C(1) + t35 * r_i_i_C(2), 0, 0, 0; 0 (t30 * r_i_i_C(1) - t27 * r_i_i_C(2) + pkin(2)) * t11 + (t10 * r_i_i_C(1) + t14 * r_i_i_C(2) - t47) * t15 (t11 * t14 + t15 * t27 + t26 * t9) * r_i_i_C(1) + (-t11 * t10 + t26 * t13 + t15 * t30) * r_i_i_C(2), 0, 0, 0;];
Ja_transl  = t4;
