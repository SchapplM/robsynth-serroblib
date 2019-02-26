% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR4
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
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:47
% EndTime: 2019-02-26 21:02:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (144->23), mult. (105->24), div. (0->0), fcn. (114->9), ass. (0->21)
t41 = r_i_i_C(3) + qJ(5);
t18 = pkin(10) + qJ(3);
t16 = qJ(4) + t18;
t13 = sin(t16);
t14 = cos(t16);
t20 = cos(pkin(11));
t28 = -r_i_i_C(1) * t20 - pkin(4);
t19 = sin(pkin(11));
t35 = r_i_i_C(2) * t19;
t24 = t41 * t13 + (-t28 - t35) * t14;
t40 = pkin(3) * cos(t18) + t24;
t39 = t13 * t35 + t41 * t14;
t37 = cos(pkin(10)) * pkin(2) + pkin(1) + t40;
t21 = sin(qJ(1));
t31 = t39 * t21;
t22 = cos(qJ(1));
t30 = t39 * t22;
t27 = t28 * t13;
t25 = t19 * r_i_i_C(1) + t20 * r_i_i_C(2) + pkin(7) + pkin(8) + qJ(2);
t23 = -pkin(3) * sin(t18) + t27;
t1 = [-t37 * t21 + t25 * t22, t21, t23 * t22 + t30, t22 * t27 + t30, t22 * t13, 0; t25 * t21 + t37 * t22, -t22, t23 * t21 + t31, t21 * t27 + t31, t21 * t13, 0; 0, 0, t40, t24, -t14, 0;];
Ja_transl  = t1;
