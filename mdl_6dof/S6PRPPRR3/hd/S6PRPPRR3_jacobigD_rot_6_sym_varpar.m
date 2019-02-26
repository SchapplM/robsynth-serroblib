% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPPRR3_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (27->18), mult. (97->48), div. (0->0), fcn. (104->10), ass. (0->23)
t178 = sin(pkin(6));
t182 = sin(qJ(5));
t193 = t178 * t182;
t181 = cos(pkin(6));
t183 = sin(qJ(2));
t192 = t181 * t183;
t185 = cos(qJ(2));
t191 = t181 * t185;
t176 = sin(pkin(11));
t179 = cos(pkin(11));
t190 = t176 * t185 - t179 * t183;
t177 = sin(pkin(10));
t180 = cos(pkin(10));
t189 = -t177 * t183 + t180 * t191;
t188 = t177 * t185 + t180 * t192;
t187 = t177 * t191 + t180 * t183;
t186 = -t177 * t192 + t180 * t185;
t184 = cos(qJ(5));
t175 = t186 * qJD(2);
t174 = t187 * qJD(2);
t173 = t188 * qJD(2);
t172 = t189 * qJD(2);
t1 = [0, 0, 0, 0, -t174 * t176 - t175 * t179 (-t174 * t179 + t175 * t176) * t182 + ((t187 * t176 + t186 * t179) * t184 - t177 * t193) * qJD(5); 0, 0, 0, 0, t172 * t176 - t173 * t179 (t172 * t179 + t173 * t176) * t182 + ((-t189 * t176 + t188 * t179) * t184 + t180 * t193) * qJD(5); 0, 0, 0, 0, t190 * t178 * qJD(2), -t181 * qJD(5) * t182 + (-t190 * qJD(5) * t184 + (t183 * t176 + t185 * t179) * t182 * qJD(2)) * t178;];
JgD_rot  = t1;
