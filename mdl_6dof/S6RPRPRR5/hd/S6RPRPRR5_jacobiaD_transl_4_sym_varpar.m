% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:16
% EndTime: 2019-02-26 20:51:16
% DurationCPUTime: 0.09s
% Computational Cost: add. (96->25), mult. (140->39), div. (0->0), fcn. (98->5), ass. (0->16)
t142 = pkin(10) + qJ(3);
t140 = sin(t142);
t141 = cos(t142);
t154 = r_i_i_C(3) + qJ(4);
t156 = pkin(3) + r_i_i_C(1);
t158 = (t156 * t140 - t154 * t141) * qJD(3) - t140 * qJD(4);
t155 = r_i_i_C(2) + pkin(7) + qJ(2);
t144 = sin(qJ(1));
t153 = qJD(1) * t144;
t145 = cos(qJ(1));
t152 = qJD(1) * t145;
t151 = qJD(3) * t141;
t149 = qJD(3) * t154;
t148 = -t156 * qJD(3) + qJD(4);
t147 = -t154 * t140 - t156 * t141 - cos(pkin(10)) * pkin(2) - pkin(1);
t1 = [t145 * qJD(2) + t158 * t144 + (-t155 * t144 + t147 * t145) * qJD(1), t152 (-t145 * t149 + t156 * t153) * t140 + (t148 * t145 - t154 * t153) * t141, -t140 * t153 + t145 * t151, 0, 0; t144 * qJD(2) - t158 * t145 + (t147 * t144 + t155 * t145) * qJD(1), t153 (-t144 * t149 - t156 * t152) * t140 + (t148 * t144 + t154 * t152) * t141, t140 * t152 + t144 * t151, 0, 0; 0, 0, -t158, qJD(3) * t140, 0, 0;];
JaD_transl  = t1;
