% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:50
% EndTime: 2019-02-26 21:00:50
% DurationCPUTime: 0.15s
% Computational Cost: add. (219->35), mult. (133->43), div. (0->0), fcn. (143->12), ass. (0->31)
t27 = qJ(3) + qJ(4);
t22 = pkin(11) + t27;
t17 = sin(t22);
t18 = cos(t22);
t28 = sin(qJ(6));
t44 = r_i_i_C(2) * t28;
t52 = pkin(9) + r_i_i_C(3);
t53 = t17 * t44 + t18 * t52;
t50 = t52 * t17 + t18 * pkin(5) + pkin(4) * cos(t27);
t29 = cos(qJ(6));
t45 = r_i_i_C(1) * t29;
t31 = (-pkin(5) - t45) * t17 - pkin(4) * sin(t27);
t24 = cos(qJ(3)) * pkin(3);
t48 = pkin(2) + t24 + t50;
t26 = qJ(1) + pkin(10);
t20 = sin(t26);
t41 = t20 * t28;
t40 = t20 * t29;
t21 = cos(t26);
t39 = t21 * t28;
t38 = t21 * t29;
t37 = t53 * t20;
t36 = t53 * t21;
t32 = -sin(qJ(3)) * pkin(3) + t31;
t30 = (-t44 + t45) * t18 + t50;
t25 = -qJ(5) - pkin(8) - pkin(7);
t4 = t18 * t38 + t41;
t3 = -t18 * t39 + t40;
t2 = -t18 * t40 + t39;
t1 = t18 * t41 + t38;
t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t21 * t25 - t48 * t20, 0, t32 * t21 + t36, t31 * t21 + t36, t20, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t20 * t25 + t48 * t21, 0, t32 * t20 + t37, t31 * t20 + t37, -t21, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 1, t24 + t30, t30, 0 (-r_i_i_C(1) * t28 - r_i_i_C(2) * t29) * t17;];
Ja_transl  = t5;
