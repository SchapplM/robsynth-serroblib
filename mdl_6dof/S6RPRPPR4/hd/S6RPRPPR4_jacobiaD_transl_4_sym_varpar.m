% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR4_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR4_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:59
% EndTime: 2019-02-26 20:40:59
% DurationCPUTime: 0.10s
% Computational Cost: add. (116->28), mult. (184->43), div. (0->0), fcn. (139->7), ass. (0->18)
t172 = pkin(9) + qJ(3);
t170 = sin(t172);
t171 = cos(t172);
t173 = sin(pkin(10));
t174 = cos(pkin(10));
t182 = r_i_i_C(1) * t174 - r_i_i_C(2) * t173 + pkin(3);
t188 = r_i_i_C(3) + qJ(4);
t190 = (t182 * t170 - t188 * t171) * qJD(3) - t170 * qJD(4);
t176 = sin(qJ(1));
t187 = qJD(1) * t176;
t177 = cos(qJ(1));
t186 = qJD(1) * t177;
t185 = qJD(3) * t171;
t183 = qJD(3) * t188;
t181 = t173 * r_i_i_C(1) + t174 * r_i_i_C(2) + pkin(7) + qJ(2);
t180 = -t182 * qJD(3) + qJD(4);
t179 = -t188 * t170 - t182 * t171 - cos(pkin(9)) * pkin(2) - pkin(1);
t1 = [t177 * qJD(2) + t190 * t176 + (-t181 * t176 + t179 * t177) * qJD(1), t186 (-t177 * t183 + t182 * t187) * t170 + (t180 * t177 - t188 * t187) * t171, -t170 * t187 + t177 * t185, 0, 0; t176 * qJD(2) - t190 * t177 + (t179 * t176 + t181 * t177) * qJD(1), t187 (-t176 * t183 - t182 * t186) * t170 + (t180 * t176 + t188 * t186) * t171, t170 * t186 + t176 * t185, 0, 0; 0, 0, -t190, qJD(3) * t170, 0, 0;];
JaD_transl  = t1;
