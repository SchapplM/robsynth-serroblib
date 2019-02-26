% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:12
% EndTime: 2019-02-26 21:03:13
% DurationCPUTime: 0.15s
% Computational Cost: add. (189->32), mult. (141->41), div. (0->0), fcn. (154->9), ass. (0->28)
t25 = sin(qJ(6));
t27 = cos(qJ(6));
t49 = r_i_i_C(1) * t25 + r_i_i_C(2) * t27;
t24 = pkin(10) + qJ(3);
t22 = qJ(4) + t24;
t19 = sin(t22);
t20 = cos(t22);
t37 = pkin(4) + pkin(9) + r_i_i_C(3);
t48 = t19 * qJ(5) + t37 * t20;
t46 = (qJ(5) + t49) * t20;
t18 = pkin(3) * cos(t24);
t45 = t18 + cos(pkin(10)) * pkin(2) + pkin(1) + t48;
t42 = pkin(5) + pkin(8) + pkin(7) + qJ(2);
t26 = sin(qJ(1));
t41 = t26 * t25;
t40 = t26 * t27;
t28 = cos(qJ(1));
t39 = t28 * t19;
t36 = t46 * t26;
t35 = t46 * t28;
t31 = t37 * t19;
t30 = t49 * t19 + t48;
t29 = -pkin(3) * sin(t24) - t31;
t4 = -t19 * t41 + t27 * t28;
t3 = t19 * t40 + t25 * t28;
t2 = t25 * t39 + t40;
t1 = t27 * t39 - t41;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t45 * t26 + t42 * t28, t26, t29 * t28 + t35, -t28 * t31 + t35, t39, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t42 * t26 + t45 * t28, -t28, t29 * t26 + t36, -t26 * t31 + t36, t26 * t19, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, t18 + t30, t30, -t20 (-r_i_i_C(1) * t27 + r_i_i_C(2) * t25) * t20;];
Ja_transl  = t5;
