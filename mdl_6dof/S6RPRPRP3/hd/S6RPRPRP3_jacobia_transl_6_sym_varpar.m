% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:51
% EndTime: 2019-02-26 20:44:51
% DurationCPUTime: 0.10s
% Computational Cost: add. (193->29), mult. (155->40), div. (0->0), fcn. (179->10), ass. (0->25)
t17 = cos(qJ(3));
t16 = sin(qJ(3));
t24 = r_i_i_C(2) + pkin(8) + qJ(4);
t19 = t24 * t16;
t7 = cos(pkin(10)) * pkin(4) + pkin(3);
t29 = t17 * t7 + pkin(2) + t19;
t12 = pkin(10) + qJ(5);
t10 = cos(t12);
t22 = r_i_i_C(3) + qJ(6);
t26 = pkin(5) + r_i_i_C(1);
t8 = sin(t12);
t28 = t26 * t10 + t22 * t8 + t7;
t13 = qJ(1) + pkin(9);
t9 = sin(t13);
t27 = t9 * t8;
t25 = t9 * t10;
t11 = cos(t13);
t23 = t11 * t17;
t20 = pkin(4) * sin(pkin(10)) + pkin(7);
t18 = -t28 * t16 + t24 * t17;
t4 = t10 * t23 + t27;
t3 = t8 * t23 - t25;
t2 = -t11 * t8 + t17 * t25;
t1 = t11 * t10 + t17 * t27;
t5 = [-sin(qJ(1)) * pkin(1) - t26 * t2 + t20 * t11 - t22 * t1 - t29 * t9, 0, t18 * t11, t11 * t16, t22 * t4 - t26 * t3, t3; cos(qJ(1)) * pkin(1) + t20 * t9 + t26 * t4 + t22 * t3 + t29 * t11, 0, t18 * t9, t9 * t16, -t26 * t1 + t22 * t2, t1; 0, 1, t28 * t17 + t19, -t17 (t22 * t10 - t26 * t8) * t16, t16 * t8;];
Ja_transl  = t5;
