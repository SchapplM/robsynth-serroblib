% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:34:43
% EndTime: 2019-02-26 20:34:43
% DurationCPUTime: 0.13s
% Computational Cost: add. (208->35), mult. (125->42), div. (0->0), fcn. (135->11), ass. (0->31)
t24 = pkin(11) + qJ(4);
t22 = qJ(5) + t24;
t17 = sin(t22);
t18 = cos(t22);
t26 = sin(qJ(6));
t41 = r_i_i_C(2) * t26;
t47 = pkin(9) + r_i_i_C(3);
t48 = t17 * t41 + t18 * t47;
t45 = t18 * pkin(5) + t47 * t17;
t16 = pkin(4) * cos(t24);
t44 = t16 + cos(pkin(11)) * pkin(3) + pkin(2) + t45;
t27 = cos(qJ(6));
t42 = r_i_i_C(1) * t27;
t25 = qJ(1) + pkin(10);
t20 = sin(t25);
t38 = t20 * t26;
t37 = t20 * t27;
t21 = cos(t25);
t36 = t21 * t26;
t35 = t21 * t27;
t34 = t48 * t20;
t32 = t48 * t21;
t30 = (-pkin(5) - t42) * t17;
t29 = (-t41 + t42) * t18 + t45;
t28 = -pkin(4) * sin(t24) + t30;
t23 = -pkin(8) - pkin(7) - qJ(3);
t4 = t18 * t35 + t38;
t3 = -t18 * t36 + t37;
t2 = -t18 * t37 + t36;
t1 = t18 * t38 + t35;
t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t21 * t23 - t44 * t20, 0, t20, t28 * t21 + t32, t21 * t30 + t32, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t20 * t23 + t44 * t21, 0, -t21, t28 * t20 + t34, t20 * t30 + t34, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 1, 0, t16 + t29, t29 (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * t17;];
Ja_transl  = t5;
