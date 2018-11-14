% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S4RPPP1_jacobiaD_transl_4_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_transl_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_transl_4_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPP1_jacobiaD_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_transl_4_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:52
% EndTime: 2018-11-14 13:45:52
% DurationCPUTime: 0.06s
% Computational Cost: add. (94->28), mult. (144->40), div. (0->0), fcn. (112->9), ass. (0->22)
t139 = sin(pkin(4));
t151 = t139 * (pkin(3) + r_i_i_C(1) + qJ(2));
t150 = -r_i_i_C(2) - qJ(3);
t141 = sin(qJ(1));
t149 = qJD(1) * t141;
t142 = cos(qJ(1));
t148 = qJD(1) * t142;
t147 = t139 * qJD(2);
t146 = -pkin(2) - r_i_i_C(3) - qJ(4);
t136 = pkin(4) + pkin(6);
t137 = pkin(4) - pkin(6);
t134 = sin(t136) / 0.2e1 - sin(t137) / 0.2e1;
t140 = cos(pkin(6));
t144 = t142 * t134 + t141 * t140;
t135 = cos(t137) / 0.2e1 + cos(t136) / 0.2e1;
t138 = sin(pkin(6));
t143 = t141 * t135 + t142 * t138;
t131 = -t134 * t149 + t140 * t148;
t130 = t143 * qJD(1);
t129 = t144 * qJD(1);
t128 = -t135 * t148 + t138 * t149;
t1 = [-t144 * qJD(4) - (-t142 * t135 + t141 * t138) * qJD(3) + t142 * t147 + t150 * t130 + t146 * t131 + (-t142 * pkin(1) - t141 * t151) * qJD(1), t139 * t148, -t128, -t129; -(t141 * t134 - t142 * t140) * qJD(4) + t143 * qJD(3) + t141 * t147 + t150 * t128 + t146 * t129 + (-t141 * pkin(1) + t142 * t151) * qJD(1), t139 * t149, t130, t131; 0, 0, 0, 0;];
JaD_transl  = t1;
