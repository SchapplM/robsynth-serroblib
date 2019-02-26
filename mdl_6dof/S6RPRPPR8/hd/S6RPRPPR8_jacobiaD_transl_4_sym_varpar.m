% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR8_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR8_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobiaD_transl_4_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:07
% EndTime: 2019-02-26 20:43:07
% DurationCPUTime: 0.08s
% Computational Cost: add. (48->21), mult. (142->35), div. (0->0), fcn. (98->4), ass. (0->17)
t134 = sin(qJ(3));
t136 = cos(qJ(3));
t146 = r_i_i_C(3) + qJ(4);
t147 = pkin(3) + r_i_i_C(1);
t151 = -(t146 * t134 + t147 * t136) * qJD(3) + t136 * qJD(4);
t150 = t147 * qJD(3) - qJD(4);
t148 = t147 * t134 - t146 * t136 + qJ(2);
t135 = sin(qJ(1));
t145 = qJD(1) * t135;
t137 = cos(qJ(1));
t144 = qJD(1) * t137;
t143 = qJD(3) * t135;
t142 = qJD(3) * t137;
t140 = qJD(1) * t147;
t139 = qJD(1) * t146;
t138 = qJD(2) + (-pkin(1) - pkin(7) - r_i_i_C(2)) * qJD(1) - t151;
t1 = [t138 * t137 - t148 * t145, t144 (t137 * t140 + t146 * t143) * t136 + (-t150 * t135 + t137 * t139) * t134, t134 * t143 - t136 * t144, 0, 0; t138 * t135 + t148 * t144, t145 (t135 * t140 - t146 * t142) * t136 + (t135 * t139 + t150 * t137) * t134, -t134 * t142 - t136 * t145, 0, 0; 0, 0, t151, qJD(3) * t136, 0, 0;];
JaD_transl  = t1;
