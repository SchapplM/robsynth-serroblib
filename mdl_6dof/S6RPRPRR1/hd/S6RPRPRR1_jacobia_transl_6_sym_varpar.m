% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:13
% EndTime: 2019-02-26 20:49:13
% DurationCPUTime: 0.12s
% Computational Cost: add. (211->35), mult. (128->43), div. (0->0), fcn. (138->12), ass. (0->31)
t25 = qJ(3) + pkin(11);
t22 = qJ(5) + t25;
t18 = sin(t22);
t19 = cos(t22);
t27 = sin(qJ(6));
t43 = r_i_i_C(2) * t27;
t49 = pkin(9) + r_i_i_C(3);
t50 = t18 * t43 + t19 * t49;
t47 = t19 * pkin(5) + t49 * t18;
t36 = pkin(4) * cos(t25) + cos(qJ(3)) * pkin(3);
t46 = pkin(2) + t36 + t47;
t28 = cos(qJ(6));
t44 = r_i_i_C(1) * t28;
t26 = qJ(1) + pkin(10);
t20 = sin(t26);
t40 = t20 * t27;
t39 = t20 * t28;
t21 = cos(t26);
t38 = t21 * t27;
t37 = t21 * t28;
t35 = t50 * t20;
t34 = t50 * t21;
t31 = (-pkin(5) - t44) * t18;
t30 = -pkin(4) * sin(t25) - sin(qJ(3)) * pkin(3) + t31;
t29 = (-t43 + t44) * t19 + t47;
t24 = -pkin(8) - qJ(4) - pkin(7);
t4 = t19 * t37 + t40;
t3 = -t19 * t38 + t39;
t2 = -t19 * t39 + t38;
t1 = t19 * t40 + t37;
t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t21 * t24 - t46 * t20, 0, t30 * t21 + t34, t20, t21 * t31 + t34, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t20 * t24 + t46 * t21, 0, t30 * t20 + t35, -t21, t20 * t31 + t35, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 1, t29 + t36, 0, t29 (-r_i_i_C(1) * t27 - r_i_i_C(2) * t28) * t18;];
Ja_transl  = t5;
