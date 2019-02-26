% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP2
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

function Ja_transl = S6RPRPRP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:15
% EndTime: 2019-02-26 20:44:15
% DurationCPUTime: 0.12s
% Computational Cost: add. (185->30), mult. (151->37), div. (0->0), fcn. (174->10), ass. (0->25)
t27 = pkin(8) + r_i_i_C(2);
t13 = qJ(3) + pkin(10);
t8 = sin(t13);
t31 = cos(qJ(3)) * pkin(3) + t27 * t8;
t10 = cos(t13);
t30 = t10 * pkin(4) + pkin(2) + t31;
t16 = sin(qJ(5));
t18 = cos(qJ(5));
t22 = r_i_i_C(3) + qJ(6);
t28 = pkin(5) + r_i_i_C(1);
t29 = t22 * t16 + t28 * t18 + pkin(4);
t14 = qJ(1) + pkin(9);
t9 = sin(t14);
t26 = t9 * t16;
t25 = t9 * t18;
t11 = cos(t14);
t24 = t11 * t16;
t23 = t11 * t18;
t19 = -sin(qJ(3)) * pkin(3) + t27 * t10 - t29 * t8;
t15 = -qJ(4) - pkin(7);
t4 = t10 * t23 + t26;
t3 = t10 * t24 - t25;
t2 = t10 * t25 - t24;
t1 = t10 * t26 + t23;
t5 = [-sin(qJ(1)) * pkin(1) - t11 * t15 - t28 * t2 - t22 * t1 - t30 * t9, 0, t19 * t11, t9, t22 * t4 - t28 * t3, t3; cos(qJ(1)) * pkin(1) - t9 * t15 + t28 * t4 + t22 * t3 + t30 * t11, 0, t19 * t9, -t11, -t28 * t1 + t22 * t2, t1; 0, 1, t29 * t10 + t31, 0 (-t28 * t16 + t22 * t18) * t8, t8 * t16;];
Ja_transl  = t5;
