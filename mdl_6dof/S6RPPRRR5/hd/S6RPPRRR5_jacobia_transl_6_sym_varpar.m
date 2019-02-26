% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:04
% EndTime: 2019-02-26 20:37:04
% DurationCPUTime: 0.12s
% Computational Cost: add. (105->31), mult. (127->39), div. (0->0), fcn. (139->8), ass. (0->29)
t22 = cos(qJ(6));
t51 = r_i_i_C(1) * t22 + pkin(5);
t18 = qJ(4) + qJ(5);
t16 = sin(t18);
t17 = cos(t18);
t49 = pkin(9) + r_i_i_C(3);
t50 = t16 * t49 + t51 * t17;
t46 = t49 * t17;
t43 = sin(qJ(4)) * pkin(4);
t45 = -pkin(5) * t16 - pkin(1) - qJ(3) - t43 + t46;
t19 = sin(qJ(6));
t40 = r_i_i_C(2) * t19;
t21 = sin(qJ(1));
t37 = t19 * t21;
t24 = cos(qJ(1));
t36 = t19 * t24;
t35 = t21 * t22;
t34 = t22 * t24;
t33 = qJ(2) - pkin(8) - pkin(7);
t31 = t17 * t40;
t30 = t50 * t21;
t29 = t50 * t24;
t28 = cos(qJ(4)) * pkin(4) - t31;
t26 = t46 + (t40 - t51) * t16;
t4 = t16 * t34 - t37;
t3 = -t16 * t36 - t35;
t2 = -t16 * t35 - t36;
t1 = t16 * t37 - t34;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t45 * t21 + t33 * t24, t21, t24, t28 * t24 + t29, -t24 * t31 + t29, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t33 * t21 - t45 * t24, -t24, t21, t28 * t21 + t30, -t21 * t31 + t30, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 0, 0, t26 - t43, t26 (-r_i_i_C(1) * t19 - r_i_i_C(2) * t22) * t17;];
Ja_transl  = t5;
