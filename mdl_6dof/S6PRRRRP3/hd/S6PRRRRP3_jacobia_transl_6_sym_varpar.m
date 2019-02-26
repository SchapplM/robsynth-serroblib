% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:20
% EndTime: 2019-02-26 20:16:20
% DurationCPUTime: 0.16s
% Computational Cost: add. (257->46), mult. (448->76), div. (0->0), fcn. (561->12), ass. (0->40)
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t28 = qJ(4) + qJ(5);
t25 = cos(t28);
t19 = pkin(5) * t25 + cos(qJ(4)) * pkin(4);
t24 = sin(t28);
t38 = r_i_i_C(1) * t25 - r_i_i_C(2) * t24 + pkin(3) + t19;
t49 = r_i_i_C(3) + qJ(6) + pkin(10) + pkin(9);
t53 = t49 * t31 + t38 * t33 + pkin(2);
t29 = sin(pkin(11));
t32 = sin(qJ(2));
t34 = cos(qJ(2));
t44 = cos(pkin(11));
t45 = cos(pkin(6));
t39 = t45 * t44;
t11 = t29 * t32 - t34 * t39;
t12 = t29 * t34 + t32 * t39;
t30 = sin(pkin(6));
t42 = t30 * t44;
t8 = t12 * t33 - t31 * t42;
t41 = t11 * t25 - t8 * t24;
t52 = t41 * r_i_i_C(1) + (-t11 * t24 - t8 * t25) * r_i_i_C(2);
t43 = t29 * t45;
t14 = -t32 * t43 + t44 * t34;
t48 = t30 * t31;
t10 = t14 * t33 + t29 * t48;
t13 = t44 * t32 + t34 * t43;
t40 = -t10 * t24 + t13 * t25;
t51 = t40 * r_i_i_C(1) + (-t10 * t25 - t13 * t24) * r_i_i_C(2);
t47 = t30 * t33;
t16 = t45 * t31 + t32 * t47;
t46 = t30 * t34;
t37 = -t16 * t24 - t25 * t46;
t50 = t37 * r_i_i_C(1) + (-t16 * t25 + t24 * t46) * r_i_i_C(2);
t18 = pkin(5) * t24 + sin(qJ(4)) * pkin(4);
t36 = t24 * r_i_i_C(1) + t25 * r_i_i_C(2) + pkin(8) + t18;
t15 = t32 * t48 - t45 * t33;
t9 = t14 * t31 - t29 * t47;
t7 = t12 * t31 + t33 * t42;
t1 = [0, -t13 * t53 + t36 * t14, t49 * t10 - t38 * t9, -t10 * t18 + t13 * t19 + t51, t40 * pkin(5) + t51, t9; 0, -t11 * t53 + t36 * t12, -t38 * t7 + t49 * t8, t11 * t19 - t8 * t18 + t52, t41 * pkin(5) + t52, t7; 1 (t36 * t32 + t53 * t34) * t30, -t38 * t15 + t49 * t16, -t16 * t18 - t19 * t46 + t50, t37 * pkin(5) + t50, t15;];
Ja_transl  = t1;
