% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
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

function JaD_transl = S4RPPP1_jacobiaD_transl_3_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_transl_3_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_transl_3_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPP1_jacobiaD_transl_3_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_transl_3_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:47
% EndTime: 2018-11-14 13:45:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (64->23), mult. (102->36), div. (0->0), fcn. (78->9), ass. (0->19)
t138 = sin(pkin(4));
t149 = t138 * (r_i_i_C(1) + qJ(2));
t147 = -r_i_i_C(3) - qJ(3);
t140 = sin(qJ(1));
t146 = qJD(1) * t140;
t141 = cos(qJ(1));
t145 = qJD(1) * t141;
t144 = t138 * qJD(2);
t143 = qJD(1) * (-pkin(2) + r_i_i_C(2));
t135 = pkin(4) + pkin(6);
t136 = pkin(4) - pkin(6);
t134 = cos(t136) / 0.2e1 + cos(t135) / 0.2e1;
t137 = sin(pkin(6));
t142 = t140 * t134 + t141 * t137;
t139 = cos(pkin(6));
t133 = sin(t135) / 0.2e1 - sin(t136) / 0.2e1;
t130 = t142 * qJD(1);
t128 = -t134 * t145 + t137 * t146;
t1 = [-(-t141 * t134 + t140 * t137) * qJD(3) + t141 * t144 + (-t133 * t140 + t139 * t141) * t143 + t147 * t130 + (-t141 * pkin(1) - t140 * t149) * qJD(1), t138 * t145, -t128, 0; t142 * qJD(3) + t140 * t144 + (t133 * t141 + t139 * t140) * t143 + t147 * t128 + (-t140 * pkin(1) + t141 * t149) * qJD(1), t138 * t146, t130, 0; 0, 0, 0, 0;];
JaD_transl  = t1;
