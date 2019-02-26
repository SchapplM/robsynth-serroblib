% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:32:53
% EndTime: 2019-02-26 22:32:53
% DurationCPUTime: 0.13s
% Computational Cost: add. (157->29), mult. (196->39), div. (0->0), fcn. (217->8), ass. (0->30)
t23 = sin(qJ(4));
t26 = cos(qJ(4));
t34 = r_i_i_C(3) + qJ(5);
t45 = pkin(4) + r_i_i_C(1);
t51 = t34 * t23 + t45 * t26;
t50 = pkin(9) + r_i_i_C(2);
t22 = qJ(2) + qJ(3);
t20 = cos(t22);
t48 = t20 * t50;
t19 = sin(t22);
t47 = t20 * pkin(3) + t50 * t19;
t21 = cos(qJ(2)) * pkin(2);
t46 = pkin(1) + t21 + t47;
t25 = sin(qJ(1));
t43 = t25 * t48;
t38 = t23 * t25;
t37 = t25 * t26;
t27 = cos(qJ(1));
t36 = t26 * t27;
t35 = t27 * t23;
t33 = t27 * t48;
t31 = t51 * t20 + t47;
t30 = (-pkin(3) - t51) * t19;
t29 = -sin(qJ(2)) * pkin(2) + t30;
t28 = -pkin(8) - pkin(7);
t4 = t20 * t36 + t38;
t3 = t20 * t35 - t37;
t2 = t20 * t37 - t35;
t1 = t20 * t38 + t36;
t5 = [-t34 * t1 - t45 * t2 - t46 * t25 - t27 * t28, t29 * t27 + t33, t27 * t30 + t33, -t45 * t3 + t34 * t4, t3, 0; -t25 * t28 + t46 * t27 + t34 * t3 + t45 * t4, t29 * t25 + t43, t25 * t30 + t43, -t45 * t1 + t34 * t2, t1, 0; 0, t21 + t31, t31 (-t45 * t23 + t34 * t26) * t19, t19 * t23, 0;];
Ja_transl  = t5;
