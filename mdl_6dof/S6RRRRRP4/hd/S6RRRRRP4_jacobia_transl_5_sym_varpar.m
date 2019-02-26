% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP4
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
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRP4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:20
% EndTime: 2019-02-26 22:41:20
% DurationCPUTime: 0.12s
% Computational Cost: add. (179->38), mult. (165->52), div. (0->0), fcn. (177->10), ass. (0->39)
t28 = cos(qJ(4));
t16 = t28 * pkin(4) + pkin(3);
t24 = qJ(2) + qJ(3);
t19 = sin(t24);
t21 = cos(t24);
t30 = -pkin(10) - pkin(9);
t54 = t21 * t16 + (r_i_i_C(3) - t30) * t19;
t22 = cos(qJ(2)) * pkin(2);
t53 = pkin(1) + t22 + t54;
t23 = qJ(4) + qJ(5);
t20 = cos(t23);
t29 = cos(qJ(1));
t40 = t29 * t20;
t18 = sin(t23);
t27 = sin(qJ(1));
t43 = t27 * t18;
t5 = t21 * t43 + t40;
t41 = t29 * t18;
t42 = t27 * t20;
t6 = -t21 * t42 + t41;
t52 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t21 * t41 + t42;
t8 = t21 * t40 + t43;
t51 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t25 = sin(qJ(4));
t50 = pkin(4) * t25;
t49 = r_i_i_C(1) * t20;
t48 = r_i_i_C(2) * t18;
t38 = t19 * t48;
t45 = t21 * t27;
t46 = r_i_i_C(3) * t45 + t27 * t38;
t44 = t21 * t29;
t39 = r_i_i_C(3) * t44 + t29 * t38;
t37 = pkin(8) + pkin(7) + t50;
t35 = -r_i_i_C(1) * t18 - r_i_i_C(2) * t20;
t34 = -t21 * t30 + (-t16 - t49) * t19;
t33 = (-t48 + t49) * t21 + t54;
t32 = -sin(qJ(2)) * pkin(2) + t34;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t53 * t27 + t37 * t29, t32 * t29 + t39, t34 * t29 + t39 (-t25 * t44 + t27 * t28) * pkin(4) + t51, t51, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t37 * t27 + t53 * t29, t32 * t27 + t46, t34 * t27 + t46 (-t25 * t45 - t28 * t29) * pkin(4) + t52, t52, 0; 0, t22 + t33, t33 (t35 - t50) * t19, t35 * t19, 0;];
Ja_transl  = t1;
