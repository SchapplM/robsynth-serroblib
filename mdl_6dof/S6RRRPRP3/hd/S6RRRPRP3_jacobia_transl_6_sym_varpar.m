% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:35
% EndTime: 2019-02-26 22:10:35
% DurationCPUTime: 0.14s
% Computational Cost: add. (244->35), mult. (210->44), div. (0->0), fcn. (234->10), ass. (0->32)
t22 = pkin(10) + qJ(5);
t17 = sin(t22);
t18 = cos(t22);
t35 = r_i_i_C(3) + qJ(6);
t44 = pkin(5) + r_i_i_C(1);
t48 = t35 * t17 + t44 * t18;
t14 = cos(pkin(10)) * pkin(4) + pkin(3);
t23 = qJ(2) + qJ(3);
t19 = sin(t23);
t20 = cos(t23);
t25 = -pkin(9) - qJ(4);
t46 = t20 * t14 + (r_i_i_C(2) - t25) * t19;
t21 = cos(qJ(2)) * pkin(2);
t45 = pkin(1) + t21 + t46;
t43 = r_i_i_C(2) * t20;
t27 = sin(qJ(1));
t39 = t27 * t17;
t38 = t27 * t18;
t28 = cos(qJ(1));
t37 = t28 * t17;
t36 = t28 * t18;
t34 = pkin(4) * sin(pkin(10)) + pkin(8) + pkin(7);
t32 = t48 * t20 + t46;
t31 = -t20 * t25 + (-t14 - t48) * t19;
t30 = -sin(qJ(2)) * pkin(2) + t31;
t11 = t28 * t43;
t10 = t27 * t43;
t4 = t20 * t36 + t39;
t3 = t20 * t37 - t38;
t2 = t20 * t38 - t37;
t1 = t20 * t39 + t36;
t5 = [-t35 * t1 - t44 * t2 - t45 * t27 + t34 * t28, t30 * t28 + t11, t31 * t28 + t11, t28 * t19, -t44 * t3 + t35 * t4, t3; t34 * t27 + t45 * t28 + t35 * t3 + t44 * t4, t30 * t27 + t10, t31 * t27 + t10, t27 * t19, -t44 * t1 + t35 * t2, t1; 0, t21 + t32, t32, -t20 (-t44 * t17 + t35 * t18) * t19, t19 * t17;];
Ja_transl  = t5;
