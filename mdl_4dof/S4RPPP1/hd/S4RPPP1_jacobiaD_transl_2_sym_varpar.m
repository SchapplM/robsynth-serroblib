% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
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

function JaD_transl = S4RPPP1_jacobiaD_transl_2_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_transl_2_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_transl_2_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPP1_jacobiaD_transl_2_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_transl_2_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:47
% EndTime: 2018-11-14 13:45:47
% DurationCPUTime: 0.04s
% Computational Cost: add. (26->17), mult. (52->29), div. (0->0), fcn. (38->9), ass. (0->13)
t117 = sin(pkin(4));
t124 = (r_i_i_C(3) + qJ(2)) * t117;
t122 = qJD(1) * t117;
t121 = t117 * qJD(2);
t120 = cos(qJ(1));
t119 = sin(qJ(1));
t118 = cos(pkin(6));
t116 = sin(pkin(6));
t115 = pkin(4) - pkin(6);
t114 = pkin(4) + pkin(6);
t113 = cos(t115) / 0.2e1 + cos(t114) / 0.2e1;
t112 = sin(t114) / 0.2e1 - sin(t115) / 0.2e1;
t1 = [t120 * t121 + ((t112 * t119 - t118 * t120) * r_i_i_C(1) + (t113 * t119 + t116 * t120) * r_i_i_C(2) - t120 * pkin(1) - t119 * t124) * qJD(1), t120 * t122, 0, 0; t119 * t121 + ((-t112 * t120 - t118 * t119) * r_i_i_C(1) + (-t113 * t120 + t116 * t119) * r_i_i_C(2) - t119 * pkin(1) + t120 * t124) * qJD(1), t119 * t122, 0, 0; 0, 0, 0, 0;];
JaD_transl  = t1;
