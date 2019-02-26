% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:25
% EndTime: 2019-02-26 21:16:26
% DurationCPUTime: 0.15s
% Computational Cost: add. (285->38), mult. (164->44), div. (0->0), fcn. (174->11), ass. (0->34)
t25 = pkin(11) + qJ(3);
t23 = qJ(4) + t25;
t22 = qJ(5) + t23;
t18 = sin(t22);
t19 = cos(t22);
t26 = sin(qJ(6));
t45 = r_i_i_C(2) * t26;
t52 = pkin(10) + r_i_i_C(3);
t53 = t18 * t45 + t19 * t52;
t50 = t19 * pkin(5) + t52 * t18;
t28 = cos(qJ(6));
t46 = r_i_i_C(1) * t28;
t34 = (-pkin(5) - t46) * t18;
t31 = t34 - pkin(4) * sin(t23);
t17 = pkin(4) * cos(t23);
t20 = pkin(3) * cos(t25);
t49 = t17 + t20 + cos(pkin(11)) * pkin(2) + pkin(1) + t50;
t29 = cos(qJ(1));
t42 = t26 * t29;
t27 = sin(qJ(1));
t41 = t27 * t26;
t40 = t27 * t28;
t39 = t29 * t28;
t38 = t53 * t27;
t36 = t53 * t29;
t33 = -pkin(3) * sin(t25) + t31;
t32 = (-t45 + t46) * t19 + t50;
t30 = t17 + t32;
t24 = -pkin(9) - pkin(8) - pkin(7) - qJ(2);
t5 = t19 * t39 + t41;
t4 = -t19 * t42 + t40;
t3 = -t19 * t40 + t42;
t2 = t19 * t41 + t39;
t1 = [t3 * r_i_i_C(1) + t2 * r_i_i_C(2) - t24 * t29 - t49 * t27, t27, t33 * t29 + t36, t31 * t29 + t36, t29 * t34 + t36, r_i_i_C(1) * t4 - r_i_i_C(2) * t5; t5 * r_i_i_C(1) + t4 * r_i_i_C(2) - t27 * t24 + t49 * t29, -t29, t33 * t27 + t38, t31 * t27 + t38, t27 * t34 + t38, -r_i_i_C(1) * t2 + r_i_i_C(2) * t3; 0, 0, t20 + t30, t30, t32 (-r_i_i_C(1) * t26 - r_i_i_C(2) * t28) * t18;];
Ja_transl  = t1;
