% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP1
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
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:41
% EndTime: 2019-02-26 20:43:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (148->32), mult. (115->40), div. (0->0), fcn. (128->10), ass. (0->25)
t26 = r_i_i_C(3) + qJ(6) + pkin(8);
t12 = qJ(3) + pkin(10);
t7 = sin(t12);
t31 = cos(qJ(3)) * pkin(3) + t26 * t7;
t18 = cos(qJ(5));
t5 = pkin(5) * t18 + pkin(4);
t9 = cos(t12);
t30 = t9 * t5 + pkin(2) + t31;
t29 = pkin(5) + r_i_i_C(1);
t13 = qJ(1) + pkin(9);
t8 = sin(t13);
t28 = t18 * t8;
t16 = sin(qJ(5));
t27 = t8 * t16;
t10 = cos(t13);
t25 = t10 * t16;
t24 = t10 * t18;
t21 = pkin(5) * t16 + pkin(7) + qJ(4);
t20 = r_i_i_C(1) * t18 - r_i_i_C(2) * t16 + t5;
t1 = t9 * t27 + t24;
t3 = -t9 * t25 + t28;
t19 = -sin(qJ(3)) * pkin(3) - t20 * t7 + t26 * t9;
t4 = t9 * t24 + t27;
t2 = -t9 * t28 + t25;
t6 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t21 * t10 - t30 * t8, 0, t19 * t10, t8, -t4 * r_i_i_C(2) + t29 * t3, t10 * t7; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t21 * t8 + t30 * t10, 0, t19 * t8, -t10, t2 * r_i_i_C(2) - t29 * t1, t8 * t7; 0, 1, t20 * t9 + t31, 0 (-r_i_i_C(2) * t18 - t29 * t16) * t7, -t9;];
Ja_transl  = t6;
