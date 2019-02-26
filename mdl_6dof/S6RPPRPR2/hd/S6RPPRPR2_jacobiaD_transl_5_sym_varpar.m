% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:17
% EndTime: 2019-02-26 20:26:17
% DurationCPUTime: 0.10s
% Computational Cost: add. (148->27), mult. (144->39), div. (0->0), fcn. (100->7), ass. (0->17)
t156 = pkin(10) + qJ(4);
t152 = sin(t156);
t154 = cos(t156);
t167 = r_i_i_C(3) + qJ(5);
t169 = pkin(4) - r_i_i_C(2);
t170 = t169 * t152 - t167 * t154;
t171 = t170 * qJD(4) - t152 * qJD(5);
t168 = r_i_i_C(1) + pkin(7) + qJ(3);
t157 = qJ(1) + pkin(9);
t153 = sin(t157);
t166 = qJD(1) * t153;
t155 = cos(t157);
t165 = qJD(1) * t155;
t164 = qJD(4) * t155;
t162 = -t167 * t152 - t169 * t154;
t160 = -cos(pkin(10)) * pkin(3) - pkin(2) + t162;
t1 = [t155 * qJD(3) + t171 * t153 + (-cos(qJ(1)) * pkin(1) - t168 * t153 + t160 * t155) * qJD(1), 0, t165 (-t167 * t164 + t169 * t166) * t152 + (-t167 * t166 + (-t169 * qJD(4) + qJD(5)) * t155) * t154, -t152 * t166 + t154 * t164, 0; t153 * qJD(3) - t171 * t155 + (-sin(qJ(1)) * pkin(1) + t168 * t155 + t160 * t153) * qJD(1), 0, t166, -t170 * t165 + (t162 * qJD(4) + qJD(5) * t154) * t153, t153 * qJD(4) * t154 + t152 * t165, 0; 0, 0, 0, -t171, qJD(4) * t152, 0;];
JaD_transl  = t1;
