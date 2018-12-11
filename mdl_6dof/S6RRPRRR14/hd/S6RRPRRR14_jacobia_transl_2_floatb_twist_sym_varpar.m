% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14_jacobia_transl_2_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_transl_2_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobia_transl_2_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_transl_2_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:20
% EndTime: 2018-12-10 18:38:20
% DurationCPUTime: 0.08s
% Computational Cost: add. (49->20), mult. (56->29), div. (0->0), fcn. (54->9), ass. (0->18)
t20 = pkin(6) - qJ(2);
t9 = pkin(6) + qJ(2);
t19 = sin(t9) / 0.2e1;
t18 = sin(pkin(6)) * (pkin(10) + r_i_i_C(3));
t17 = sin(t20);
t12 = sin(qJ(1));
t13 = cos(qJ(2));
t14 = cos(qJ(1));
t3 = t19 - t17 / 0.2e1;
t16 = t12 * t3 - t14 * t13;
t15 = -t12 * t13 - t14 * t3;
t11 = sin(qJ(2));
t8 = cos(t20);
t7 = cos(t9) / 0.2e1;
t4 = t8 / 0.2e1 + t7;
t2 = -t14 * t11 - t12 * t4;
t1 = t12 * t11 - t14 * t4;
t5 = [-t12 * pkin(1) + t15 * r_i_i_C(1) + t1 * r_i_i_C(2) + t14 * t18, t2 * r_i_i_C(1) + r_i_i_C(2) * t16, 0, 0, 0, 0; t14 * pkin(1) - t16 * r_i_i_C(1) + t2 * r_i_i_C(2) + t12 * t18, -t1 * r_i_i_C(1) + r_i_i_C(2) * t15, 0, 0, 0, 0; 0 (t19 + t17 / 0.2e1) * r_i_i_C(1) + (t7 - t8 / 0.2e1) * r_i_i_C(2), 0, 0, 0, 0;];
Ja_transl  = t5;
