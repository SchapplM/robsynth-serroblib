% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:42:33
% EndTime: 2019-02-26 20:42:34
% DurationCPUTime: 0.09s
% Computational Cost: add. (108->23), mult. (172->29), div. (0->0), fcn. (119->6), ass. (0->18)
t148 = qJ(3) + pkin(9);
t145 = sin(t148);
t146 = cos(t148);
t163 = r_i_i_C(3) + qJ(5);
t166 = pkin(4) - r_i_i_C(2);
t157 = t166 * t145 - t163 * t146 + sin(qJ(3)) * pkin(3);
t171 = qJD(4) + (qJ(2) + t157) * qJD(1);
t169 = t157 * qJD(3) - qJD(5) * t145;
t156 = t163 * t145 + t166 * t146 + cos(qJ(3)) * pkin(3);
t168 = -t156 * qJD(3) + t146 * qJD(5);
t151 = sin(qJ(1));
t162 = qJD(1) * t151;
t153 = cos(qJ(1));
t147 = qJD(1) * t153;
t161 = qJD(3) * t145;
t155 = qJD(1) * t156;
t154 = qJD(2) + (-pkin(1) - r_i_i_C(1) - qJ(4) - pkin(7)) * qJD(1) - t168;
t1 = [-t171 * t151 + t154 * t153, t147, -t169 * t151 + t153 * t155, -t162, -t146 * t147 + t151 * t161, 0; t154 * t151 + t171 * t153, t162, t151 * t155 + t169 * t153, t147, -t146 * t162 - t153 * t161, 0; 0, 0, t168, 0, qJD(3) * t146, 0;];
JaD_transl  = t1;
