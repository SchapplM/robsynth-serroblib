% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR11_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:51
% EndTime: 2019-02-26 22:21:51
% DurationCPUTime: 0.05s
% Computational Cost: add. (19->17), mult. (52->30), div. (0->0), fcn. (81->10), ass. (0->23)
t177 = sin(pkin(6));
t181 = sin(qJ(2));
t194 = t177 * t181;
t185 = cos(qJ(2));
t193 = t177 * t185;
t182 = sin(qJ(1));
t192 = t182 * t177;
t191 = t182 * t181;
t190 = t182 * t185;
t186 = cos(qJ(1));
t189 = t186 * t177;
t188 = t186 * t181;
t187 = t186 * t185;
t184 = cos(qJ(3));
t183 = cos(qJ(5));
t180 = sin(qJ(3));
t179 = sin(qJ(5));
t178 = cos(pkin(6));
t176 = -t178 * t191 + t187;
t175 = t178 * t190 + t188;
t174 = t178 * t188 + t190;
t173 = -t178 * t187 + t191;
t1 = [0, t192, t175, 0, -t175 (t176 * t184 + t180 * t192) * t179 - (t176 * t180 - t184 * t192) * t183; 0, -t189, t173, 0, -t173 (t174 * t184 - t180 * t189) * t179 - (t174 * t180 + t184 * t189) * t183; 1, t178, -t193, 0, t193 (t178 * t180 + t184 * t194) * t179 - (-t178 * t184 + t180 * t194) * t183;];
Jg_rot  = t1;
