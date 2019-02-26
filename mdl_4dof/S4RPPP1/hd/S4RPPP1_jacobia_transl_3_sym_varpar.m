% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
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

function Ja_transl = S4RPPP1_jacobia_transl_3_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobia_transl_3_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPP1_jacobia_transl_3_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobia_transl_3_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:47
% EndTime: 2018-11-14 13:45:47
% DurationCPUTime: 0.08s
% Computational Cost: add. (56->19), mult. (62->24), div. (0->0), fcn. (67->10), ass. (0->17)
t18 = pkin(2) - r_i_i_C(2);
t17 = r_i_i_C(3) + qJ(3);
t12 = sin(pkin(4));
t16 = t12 * (r_i_i_C(1) + qJ(2));
t15 = cos(qJ(1));
t14 = sin(qJ(1));
t13 = cos(pkin(6));
t11 = sin(pkin(6));
t10 = pkin(4) - pkin(6);
t9 = pkin(4) + pkin(6);
t8 = sin(t9);
t7 = -sin(t10) / 0.2e1;
t6 = cos(t10) / 0.2e1 + cos(t9) / 0.2e1;
t5 = t8 / 0.2e1 + t7;
t3 = t15 * t11 + t14 * t6;
t1 = t14 * t11 - t15 * t6;
t2 = [-t14 * pkin(1) - t18 * (t14 * t13 + t15 * t5) + t15 * t16 - t17 * t1, t14 * t12, t3, 0; t15 * pkin(1) + t18 * (t15 * t13 - t14 * t5) + t17 * t3 + t14 * t16, -t15 * t12, t1, 0; 0, cos(pkin(4)) -t8 / 0.2e1 + t7, 0;];
Ja_transl  = t2;
