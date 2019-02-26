% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S4RPPP1_jacobia_transl_4_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobia_transl_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPP1_jacobia_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobia_transl_4_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:52
% EndTime: 2018-11-14 13:45:52
% DurationCPUTime: 0.07s
% Computational Cost: add. (79->21), mult. (86->25), div. (0->0), fcn. (93->10), ass. (0->21)
t20 = r_i_i_C(2) + qJ(3);
t19 = pkin(2) + r_i_i_C(3) + qJ(4);
t14 = sin(pkin(4));
t18 = t14 * (pkin(3) + r_i_i_C(1) + qJ(2));
t17 = cos(qJ(1));
t16 = sin(qJ(1));
t15 = cos(pkin(6));
t13 = sin(pkin(6));
t12 = pkin(4) - pkin(6);
t11 = pkin(4) + pkin(6);
t10 = cos(t11);
t9 = sin(t11);
t8 = cos(t12) / 0.2e1;
t7 = -sin(t12) / 0.2e1;
t6 = t8 + t10 / 0.2e1;
t5 = t9 / 0.2e1 + t7;
t4 = t15 * t17 - t16 * t5;
t3 = t13 * t17 + t16 * t6;
t2 = t16 * t15 + t17 * t5;
t1 = t16 * t13 - t17 * t6;
t21 = [-t16 * pkin(1) - t20 * t1 + t17 * t18 - t19 * t2, t16 * t14, t3, t4; t17 * pkin(1) + t16 * t18 + t19 * t4 + t20 * t3, -t17 * t14, t1, t2; 0, cos(pkin(4)) -t9 / 0.2e1 + t7, t8 - t10 / 0.2e1;];
Ja_transl  = t21;
