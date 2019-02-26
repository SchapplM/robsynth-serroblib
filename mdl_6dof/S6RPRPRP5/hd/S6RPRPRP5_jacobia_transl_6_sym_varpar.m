% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP5
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
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:52
% EndTime: 2019-02-26 20:45:52
% DurationCPUTime: 0.11s
% Computational Cost: add. (185->30), mult. (155->39), div. (0->0), fcn. (181->9), ass. (0->25)
t14 = pkin(9) + qJ(3);
t12 = cos(t14);
t10 = sin(t14);
t27 = r_i_i_C(2) + pkin(8) + qJ(4);
t22 = t27 * t10;
t7 = cos(pkin(10)) * pkin(4) + pkin(3);
t31 = t22 + t12 * t7 + cos(pkin(9)) * pkin(2) + pkin(1);
t13 = pkin(10) + qJ(5);
t11 = cos(t13);
t24 = r_i_i_C(3) + qJ(6);
t29 = pkin(5) + r_i_i_C(1);
t9 = sin(t13);
t30 = t29 * t11 + t24 * t9 + t7;
t18 = sin(qJ(1));
t28 = t18 * t9;
t19 = cos(qJ(1));
t26 = t12 * t19;
t25 = t18 * t11;
t21 = pkin(4) * sin(pkin(10)) + pkin(7) + qJ(2);
t20 = -t30 * t10 + t27 * t12;
t4 = t11 * t26 + t28;
t3 = t9 * t26 - t25;
t2 = t12 * t25 - t19 * t9;
t1 = t11 * t19 + t12 * t28;
t5 = [-t24 * t1 - t31 * t18 + t21 * t19 - t29 * t2, t18, t20 * t19, t19 * t10, t24 * t4 - t29 * t3, t3; t21 * t18 + t31 * t19 + t24 * t3 + t29 * t4, -t19, t20 * t18, t18 * t10, -t29 * t1 + t24 * t2, t1; 0, 0, t30 * t12 + t22, -t12 (t24 * t11 - t29 * t9) * t10, t10 * t9;];
Ja_transl  = t5;
