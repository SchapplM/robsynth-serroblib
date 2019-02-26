% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:56
% EndTime: 2019-02-26 20:46:56
% DurationCPUTime: 0.10s
% Computational Cost: add. (100->31), mult. (117->38), div. (0->0), fcn. (132->8), ass. (0->24)
t29 = pkin(5) + r_i_i_C(1);
t26 = r_i_i_C(3) + qJ(6) + pkin(8);
t8 = qJ(3) + pkin(9);
t7 = cos(t8);
t28 = t26 * t7 - sin(qJ(3)) * pkin(3);
t11 = sin(qJ(5));
t14 = cos(qJ(5));
t5 = pkin(5) * t14 + pkin(4);
t18 = r_i_i_C(1) * t14 - r_i_i_C(2) * t11 + t5;
t6 = sin(t8);
t27 = t18 * t7 + t26 * t6 + cos(qJ(3)) * pkin(3);
t16 = cos(qJ(1));
t23 = t11 * t16;
t13 = sin(qJ(1));
t22 = t13 * t11;
t21 = t13 * t14;
t20 = t14 * t16;
t19 = pkin(5) * t11 + pkin(1) + pkin(7) + qJ(4);
t3 = t6 * t23 + t21;
t1 = -t6 * t22 + t20;
t17 = t6 * t5 + qJ(2) - t28;
t4 = t6 * t20 - t22;
t2 = t6 * t21 + t23;
t9 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t19 * t13 + t17 * t16, t13, t27 * t13, t16, -t2 * r_i_i_C(2) + t29 * t1, -t13 * t7; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t17 * t13 + t19 * t16, -t16, -t27 * t16, t13, t4 * r_i_i_C(2) + t29 * t3, t16 * t7; 0, 0, -t18 * t6 + t28, 0 (-r_i_i_C(2) * t14 - t29 * t11) * t7, t6;];
Ja_transl  = t9;
