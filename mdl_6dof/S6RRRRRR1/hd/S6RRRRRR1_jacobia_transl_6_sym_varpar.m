% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:47:03
% EndTime: 2019-02-26 22:47:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (359->38), mult. (206->47), div. (0->0), fcn. (214->12), ass. (0->37)
t28 = qJ(2) + qJ(3);
t25 = qJ(4) + t28;
t23 = qJ(5) + t25;
t19 = sin(t23);
t20 = cos(t23);
t29 = sin(qJ(6));
t50 = r_i_i_C(2) * t29;
t57 = pkin(11) + r_i_i_C(3);
t58 = t19 * t50 + t20 * t57;
t31 = cos(qJ(6));
t51 = r_i_i_C(1) * t31;
t39 = (-pkin(5) - t51) * t19;
t35 = t39 - pkin(4) * sin(t25);
t55 = t20 * pkin(5) + t57 * t19;
t37 = -pkin(3) * sin(t28) + t35;
t18 = pkin(4) * cos(t25);
t21 = pkin(3) * cos(t28);
t27 = cos(qJ(2)) * pkin(2);
t54 = pkin(1) + t18 + t21 + t27 + t55;
t30 = sin(qJ(1));
t47 = t30 * t29;
t46 = t30 * t31;
t32 = cos(qJ(1));
t45 = t32 * t29;
t44 = t32 * t31;
t42 = t58 * t30;
t41 = t58 * t32;
t38 = -sin(qJ(2)) * pkin(2) + t37;
t36 = (-t50 + t51) * t20 + t55;
t34 = t18 + t36;
t33 = t21 + t34;
t26 = -pkin(10) - pkin(9) - pkin(8) - pkin(7);
t4 = t20 * t44 + t47;
t3 = -t20 * t45 + t46;
t2 = -t20 * t46 + t45;
t1 = t20 * t47 + t44;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t32 * t26 - t30 * t54, t38 * t32 + t41, t37 * t32 + t41, t35 * t32 + t41, t32 * t39 + t41, t3 * r_i_i_C(1) - t4 * r_i_i_C(2); t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t30 * t26 + t32 * t54, t38 * t30 + t42, t37 * t30 + t42, t35 * t30 + t42, t30 * t39 + t42, -t1 * r_i_i_C(1) + r_i_i_C(2) * t2; 0, t27 + t33, t33, t34, t36 (-r_i_i_C(1) * t29 - r_i_i_C(2) * t31) * t19;];
Ja_transl  = t5;
