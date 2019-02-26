% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function Ja_transl = S6RRRRPP2_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP2_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobia_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:25:31
% EndTime: 2019-02-26 22:25:31
% DurationCPUTime: 0.12s
% Computational Cost: add. (100->28), mult. (121->39), div. (0->0), fcn. (129->8), ass. (0->29)
t19 = qJ(2) + qJ(3);
t16 = sin(t19);
t17 = cos(t19);
t20 = sin(qJ(4));
t39 = r_i_i_C(2) * t20;
t45 = pkin(9) + r_i_i_C(3);
t46 = t16 * t39 + t17 * t45;
t43 = t17 * pkin(3) + t45 * t16;
t18 = cos(qJ(2)) * pkin(2);
t42 = pkin(1) + t18 + t43;
t23 = cos(qJ(4));
t40 = r_i_i_C(1) * t23;
t22 = sin(qJ(1));
t36 = t20 * t22;
t24 = cos(qJ(1));
t35 = t20 * t24;
t34 = t22 * t23;
t33 = t23 * t24;
t32 = t46 * t22;
t30 = t46 * t24;
t28 = (-pkin(3) - t40) * t16;
t27 = (-t39 + t40) * t17 + t43;
t26 = -sin(qJ(2)) * pkin(2) + t28;
t25 = -pkin(8) - pkin(7);
t4 = t17 * t33 + t36;
t3 = -t17 * t35 + t34;
t2 = -t17 * t34 + t35;
t1 = t17 * t36 + t33;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t42 * t22 - t24 * t25, t26 * t24 + t30, t24 * t28 + t30, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t22 * t25 + t42 * t24, t26 * t22 + t32, t22 * t28 + t32, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t18 + t27, t27 (-r_i_i_C(1) * t20 - r_i_i_C(2) * t23) * t16, 0, 0;];
Ja_transl  = t5;
