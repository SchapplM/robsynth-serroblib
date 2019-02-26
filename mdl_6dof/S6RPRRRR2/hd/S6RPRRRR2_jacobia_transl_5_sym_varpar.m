% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR2
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

function Ja_transl = S6RPRRRR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:15:25
% EndTime: 2019-02-26 21:15:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (150->31), mult. (123->44), div. (0->0), fcn. (131->10), ass. (0->26)
t22 = qJ(3) + qJ(4);
t18 = sin(t22);
t19 = cos(t22);
t23 = sin(qJ(5));
t42 = pkin(9) + r_i_i_C(3);
t43 = r_i_i_C(2) * t18 * t23 + t19 * t42;
t40 = t19 * pkin(4) + t42 * t18;
t20 = cos(qJ(3)) * pkin(3);
t39 = pkin(2) + t20 + t40;
t35 = t19 * t23;
t25 = cos(qJ(5));
t34 = t19 * t25;
t21 = qJ(1) + pkin(11);
t16 = sin(t21);
t33 = t43 * t16;
t17 = cos(t21);
t32 = t43 * t17;
t29 = (-r_i_i_C(1) * t25 - pkin(4)) * t18;
t28 = r_i_i_C(1) * t34 - r_i_i_C(2) * t35 + t40;
t27 = -sin(qJ(3)) * pkin(3) + t29;
t26 = -pkin(8) - pkin(7);
t4 = t16 * t23 + t17 * t34;
t3 = t16 * t25 - t17 * t35;
t2 = -t16 * t34 + t17 * t23;
t1 = t16 * t35 + t17 * t25;
t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t17 * t26 - t39 * t16, 0, t17 * t27 + t32, t17 * t29 + t32, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t16 * t26 + t39 * t17, 0, t16 * t27 + t33, t16 * t29 + t33, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 1, t20 + t28, t28 (-r_i_i_C(1) * t23 - r_i_i_C(2) * t25) * t18, 0;];
Ja_transl  = t5;
