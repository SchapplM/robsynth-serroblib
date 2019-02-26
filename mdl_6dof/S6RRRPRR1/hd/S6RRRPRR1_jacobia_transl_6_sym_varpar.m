% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:03
% EndTime: 2019-02-26 22:16:03
% DurationCPUTime: 0.15s
% Computational Cost: add. (294->37), mult. (170->45), div. (0->0), fcn. (180->12), ass. (0->34)
t28 = qJ(2) + qJ(3);
t24 = pkin(11) + t28;
t23 = qJ(5) + t24;
t19 = sin(t23);
t20 = cos(t23);
t29 = sin(qJ(6));
t49 = r_i_i_C(2) * t29;
t55 = pkin(10) + r_i_i_C(3);
t56 = t19 * t49 + t20 * t55;
t53 = t20 * pkin(5) + t55 * t19;
t31 = cos(qJ(6));
t50 = r_i_i_C(1) * t31;
t37 = (-pkin(5) - t50) * t19;
t35 = -pkin(3) * sin(t28) - pkin(4) * sin(t24) + t37;
t27 = cos(qJ(2)) * pkin(2);
t42 = pkin(4) * cos(t24) + pkin(3) * cos(t28);
t52 = pkin(1) + t27 + t42 + t53;
t30 = sin(qJ(1));
t46 = t30 * t29;
t45 = t30 * t31;
t32 = cos(qJ(1));
t44 = t32 * t29;
t43 = t32 * t31;
t40 = t56 * t30;
t39 = t56 * t32;
t36 = -sin(qJ(2)) * pkin(2) + t35;
t34 = (-t49 + t50) * t20 + t53;
t33 = t34 + t42;
t26 = -pkin(9) - qJ(4) - pkin(8) - pkin(7);
t4 = t20 * t43 + t46;
t3 = -t20 * t44 + t45;
t2 = -t20 * t45 + t44;
t1 = t20 * t46 + t43;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t32 * t26 - t52 * t30, t32 * t36 + t39, t32 * t35 + t39, t30, t32 * t37 + t39, t3 * r_i_i_C(1) - t4 * r_i_i_C(2); t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t30 * t26 + t52 * t32, t30 * t36 + t40, t30 * t35 + t40, -t32, t30 * t37 + t40, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2); 0, t27 + t33, t33, 0, t34 (-r_i_i_C(1) * t29 - r_i_i_C(2) * t31) * t19;];
Ja_transl  = t5;
