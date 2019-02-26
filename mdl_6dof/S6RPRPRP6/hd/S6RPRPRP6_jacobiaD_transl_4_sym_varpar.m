% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP6_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP6_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:32
% EndTime: 2019-02-26 20:46:32
% DurationCPUTime: 0.08s
% Computational Cost: add. (96->25), mult. (140->39), div. (0->0), fcn. (98->5), ass. (0->16)
t146 = pkin(9) + qJ(3);
t144 = sin(t146);
t145 = cos(t146);
t158 = r_i_i_C(3) + qJ(4);
t160 = pkin(3) - r_i_i_C(2);
t162 = (t160 * t144 - t158 * t145) * qJD(3) - t144 * qJD(4);
t159 = r_i_i_C(1) + pkin(7) + qJ(2);
t148 = sin(qJ(1));
t157 = qJD(1) * t148;
t149 = cos(qJ(1));
t156 = qJD(1) * t149;
t155 = qJD(3) * t145;
t153 = qJD(3) * t158;
t152 = -t160 * qJD(3) + qJD(4);
t151 = -t158 * t144 - t160 * t145 - cos(pkin(9)) * pkin(2) - pkin(1);
t1 = [t149 * qJD(2) + t162 * t148 + (-t159 * t148 + t151 * t149) * qJD(1), t156 (-t149 * t153 + t160 * t157) * t144 + (t152 * t149 - t158 * t157) * t145, -t144 * t157 + t149 * t155, 0, 0; t148 * qJD(2) - t162 * t149 + (t151 * t148 + t159 * t149) * qJD(1), t157 (-t148 * t153 - t160 * t156) * t144 + (t152 * t148 + t158 * t156) * t145, t144 * t156 + t148 * t155, 0, 0; 0, 0, -t162, qJD(3) * t144, 0, 0;];
JaD_transl  = t1;
