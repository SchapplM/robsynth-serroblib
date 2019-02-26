% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:25:30
% EndTime: 2019-02-26 22:25:31
% DurationCPUTime: 0.14s
% Computational Cost: add. (204->33), mult. (246->42), div. (0->0), fcn. (274->8), ass. (0->30)
t21 = sin(qJ(4));
t24 = cos(qJ(4));
t32 = pkin(4) + pkin(5) + r_i_i_C(1);
t34 = r_i_i_C(2) + qJ(5);
t46 = t34 * t21 + t32 * t24;
t20 = qJ(2) + qJ(3);
t17 = sin(t20);
t18 = cos(t20);
t33 = -r_i_i_C(3) - qJ(6);
t43 = t18 * pkin(3) + (pkin(9) + t33) * t17;
t19 = cos(qJ(2)) * pkin(2);
t42 = pkin(1) + t19 + t43;
t41 = pkin(9) * t18;
t23 = sin(qJ(1));
t38 = t23 * t21;
t37 = t23 * t24;
t25 = cos(qJ(1));
t36 = t25 * t21;
t35 = t25 * t24;
t29 = t46 * t18 + t43;
t28 = t33 * t18 + (-pkin(3) - t46) * t17;
t27 = -sin(qJ(2)) * pkin(2) + t28;
t26 = -pkin(8) - pkin(7);
t11 = t25 * t41;
t6 = t23 * t41;
t4 = t18 * t35 + t38;
t3 = t18 * t36 - t37;
t2 = t18 * t37 - t36;
t1 = t18 * t38 + t35;
t5 = [-t34 * t1 - t32 * t2 - t42 * t23 - t25 * t26, t27 * t25 + t11, t28 * t25 + t11, -t32 * t3 + t34 * t4, t3, -t25 * t17; -t23 * t26 + t42 * t25 + t34 * t3 + t32 * t4, t27 * t23 + t6, t28 * t23 + t6, -t32 * t1 + t34 * t2, t1, -t23 * t17; 0, t19 + t29, t29 (-t32 * t21 + t34 * t24) * t17, t17 * t21, t18;];
Ja_transl  = t5;
