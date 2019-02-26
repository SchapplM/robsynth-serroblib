% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:41
% EndTime: 2019-02-26 22:40:41
% DurationCPUTime: 0.16s
% Computational Cost: add. (234->45), mult. (194->57), div. (0->0), fcn. (209->10), ass. (0->38)
t28 = qJ(4) + qJ(5);
t23 = cos(t28);
t14 = pkin(5) * t23 + cos(qJ(4)) * pkin(4);
t12 = pkin(3) + t14;
t29 = qJ(2) + qJ(3);
t22 = sin(t29);
t24 = cos(t29);
t27 = -qJ(6) - pkin(10) - pkin(9);
t55 = t24 * t12 + (r_i_i_C(3) - t27) * t22;
t26 = cos(qJ(2)) * pkin(2);
t54 = pkin(1) + t26 + t55;
t32 = cos(qJ(1));
t42 = t32 * t23;
t21 = sin(t28);
t31 = sin(qJ(1));
t45 = t31 * t21;
t5 = t24 * t45 + t42;
t43 = t32 * t21;
t44 = t31 * t23;
t6 = -t24 * t44 + t43;
t53 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t24 * t43 + t44;
t8 = t24 * t42 + t45;
t52 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t51 = r_i_i_C(1) * t23;
t50 = r_i_i_C(2) * t21;
t49 = r_i_i_C(2) * t23;
t47 = t24 * t31;
t46 = t24 * t32;
t38 = t22 * t50;
t41 = r_i_i_C(3) * t47 + t31 * t38;
t40 = r_i_i_C(3) * t46 + t32 * t38;
t13 = pkin(5) * t21 + sin(qJ(4)) * pkin(4);
t39 = t13 + pkin(8) + pkin(7);
t36 = -t24 * t27 + (-t12 - t51) * t22;
t35 = (-t50 + t51) * t24 + t55;
t34 = -sin(qJ(2)) * pkin(2) + t36;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t54 * t31 + t39 * t32, t32 * t34 + t40, t32 * t36 + t40, -t13 * t46 + t31 * t14 + t52, t7 * pkin(5) + t52, t32 * t22; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t39 * t31 + t54 * t32, t31 * t34 + t41, t31 * t36 + t41, -t13 * t47 - t32 * t14 + t53, -pkin(5) * t5 + t53, t31 * t22; 0, t26 + t35, t35 (-r_i_i_C(1) * t21 - t13 - t49) * t22 (-t49 + (-pkin(5) - r_i_i_C(1)) * t21) * t22, -t24;];
Ja_transl  = t1;
