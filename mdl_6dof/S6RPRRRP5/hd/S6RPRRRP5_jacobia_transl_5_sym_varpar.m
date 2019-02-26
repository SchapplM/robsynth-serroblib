% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:23
% EndTime: 2019-02-26 21:10:24
% DurationCPUTime: 0.11s
% Computational Cost: add. (156->32), mult. (123->40), div. (0->0), fcn. (133->9), ass. (0->30)
t22 = pkin(10) + qJ(3);
t20 = qJ(4) + t22;
t17 = sin(t20);
t18 = cos(t20);
t23 = sin(qJ(5));
t40 = r_i_i_C(2) * t23;
t46 = pkin(9) + r_i_i_C(3);
t47 = t17 * t40 + t18 * t46;
t44 = t18 * pkin(4) + t46 * t17;
t16 = pkin(3) * cos(t22);
t43 = t16 + cos(pkin(10)) * pkin(2) + pkin(1) + t44;
t25 = cos(qJ(5));
t41 = r_i_i_C(1) * t25;
t26 = cos(qJ(1));
t37 = t23 * t26;
t24 = sin(qJ(1));
t36 = t24 * t23;
t35 = t24 * t25;
t34 = t25 * t26;
t33 = t47 * t24;
t31 = t47 * t26;
t29 = (-pkin(4) - t41) * t17;
t28 = (-t40 + t41) * t18 + t44;
t27 = -pkin(3) * sin(t22) + t29;
t21 = -pkin(8) - pkin(7) - qJ(2);
t4 = t18 * t34 + t36;
t3 = -t18 * t37 + t35;
t2 = -t18 * t35 + t37;
t1 = t18 * t36 + t34;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t21 * t26 - t43 * t24, t24, t27 * t26 + t31, t26 * t29 + t31, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t24 * t21 + t43 * t26, -t26, t27 * t24 + t33, t24 * t29 + t33, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 0, t16 + t28, t28 (-r_i_i_C(1) * t23 - r_i_i_C(2) * t25) * t17, 0;];
Ja_transl  = t5;
