% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPP8_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:35
% EndTime: 2019-02-26 22:29:35
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
t175 = sin(pkin(6));
t180 = cos(qJ(3));
t194 = t175 * t180;
t182 = cos(qJ(1));
t193 = t175 * t182;
t178 = sin(qJ(2));
t179 = sin(qJ(1));
t192 = t178 * t179;
t191 = t178 * t182;
t181 = cos(qJ(2));
t190 = t179 * t181;
t189 = t182 * t181;
t188 = qJD(1) * t175;
t177 = sin(qJ(3));
t187 = qJD(2) * t177;
t176 = cos(pkin(6));
t186 = t176 * t189 - t192;
t185 = t176 * t190 + t191;
t184 = t176 * t191 + t190;
t183 = -t176 * t192 + t189;
t1 = [0, t182 * t188, t186 * qJD(1) + t183 * qJD(2) (t179 * t175 * t177 + t183 * t180) * qJD(3) - t185 * t187 + (-t184 * t177 - t180 * t193) * qJD(1), 0, 0; 0, t179 * t188, t185 * qJD(1) + t184 * qJD(2) (-t177 * t193 + t184 * t180) * qJD(3) + t186 * t187 + (t183 * t177 - t179 * t194) * qJD(1), 0, 0; 0, 0, t175 * qJD(2) * t178, t175 * t181 * t187 + (t176 * t177 + t178 * t194) * qJD(3), 0, 0;];
JgD_rot  = t1;
