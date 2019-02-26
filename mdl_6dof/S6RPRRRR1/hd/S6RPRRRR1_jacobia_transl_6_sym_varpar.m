% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR1
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
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:50
% EndTime: 2019-02-26 21:14:51
% DurationCPUTime: 0.15s
% Computational Cost: add. (268->36), mult. (164->45), div. (0->0), fcn. (172->12), ass. (0->34)
t27 = qJ(3) + qJ(4);
t23 = qJ(5) + t27;
t18 = sin(t23);
t19 = cos(t23);
t28 = sin(qJ(6));
t45 = r_i_i_C(2) * t28;
t52 = pkin(10) + r_i_i_C(3);
t53 = t18 * t45 + t19 * t52;
t50 = t19 * pkin(5) + t52 * t18;
t29 = cos(qJ(6));
t46 = r_i_i_C(1) * t29;
t34 = (-pkin(5) - t46) * t18;
t31 = t34 - pkin(4) * sin(t27);
t17 = pkin(4) * cos(t27);
t24 = cos(qJ(3)) * pkin(3);
t49 = pkin(2) + t17 + t24 + t50;
t25 = qJ(1) + pkin(11);
t20 = sin(t25);
t42 = t20 * t28;
t41 = t20 * t29;
t21 = cos(t25);
t40 = t21 * t28;
t39 = t21 * t29;
t38 = t53 * t20;
t37 = t53 * t21;
t33 = -sin(qJ(3)) * pkin(3) + t31;
t32 = (-t45 + t46) * t19 + t50;
t30 = t17 + t32;
t26 = -pkin(9) - pkin(8) - pkin(7);
t4 = t19 * t39 + t42;
t3 = -t19 * t40 + t41;
t2 = -t19 * t41 + t40;
t1 = t19 * t42 + t39;
t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t21 * t26 - t49 * t20, 0, t33 * t21 + t37, t31 * t21 + t37, t21 * t34 + t37, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t20 * t26 + t49 * t21, 0, t33 * t20 + t38, t31 * t20 + t38, t20 * t34 + t38, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 1, t24 + t30, t30, t32 (-r_i_i_C(1) * t28 - r_i_i_C(2) * t29) * t18;];
Ja_transl  = t5;
