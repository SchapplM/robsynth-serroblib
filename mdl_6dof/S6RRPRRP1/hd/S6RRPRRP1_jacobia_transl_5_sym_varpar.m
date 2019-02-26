% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:46:00
% EndTime: 2019-02-26 21:46:00
% DurationCPUTime: 0.12s
% Computational Cost: add. (159->32), mult. (126->41), div. (0->0), fcn. (136->10), ass. (0->30)
t23 = qJ(2) + pkin(10);
t20 = qJ(4) + t23;
t18 = sin(t20);
t19 = cos(t20);
t24 = sin(qJ(5));
t42 = r_i_i_C(2) * t24;
t48 = pkin(9) + r_i_i_C(3);
t49 = t18 * t42 + t19 * t48;
t46 = t19 * pkin(4) + t48 * t18;
t35 = pkin(3) * cos(t23) + cos(qJ(2)) * pkin(2);
t45 = pkin(1) + t35 + t46;
t26 = cos(qJ(5));
t43 = r_i_i_C(1) * t26;
t27 = cos(qJ(1));
t39 = t24 * t27;
t25 = sin(qJ(1));
t38 = t25 * t24;
t37 = t25 * t26;
t36 = t26 * t27;
t34 = t49 * t25;
t32 = t49 * t27;
t30 = (-pkin(4) - t43) * t18;
t29 = -pkin(3) * sin(t23) - sin(qJ(2)) * pkin(2) + t30;
t28 = (-t42 + t43) * t19 + t46;
t22 = -pkin(8) - qJ(3) - pkin(7);
t4 = t19 * t36 + t38;
t3 = -t19 * t39 + t37;
t2 = -t19 * t37 + t39;
t1 = t19 * t38 + t36;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t22 * t27 - t45 * t25, t29 * t27 + t32, t25, t27 * t30 + t32, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t25 * t22 + t45 * t27, t29 * t25 + t34, -t27, t25 * t30 + t34, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t28 + t35, 0, t28 (-r_i_i_C(1) * t24 - r_i_i_C(2) * t26) * t18, 0;];
Ja_transl  = t5;
