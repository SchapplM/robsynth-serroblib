% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:52
% EndTime: 2019-02-26 20:19:52
% DurationCPUTime: 0.09s
% Computational Cost: add. (19->9), mult. (54->20), div. (0->0), fcn. (86->8), ass. (0->16)
t107 = sin(pkin(12));
t108 = sin(pkin(6));
t118 = t107 * t108;
t109 = cos(pkin(12));
t117 = t109 * t108;
t110 = cos(pkin(6));
t112 = sin(qJ(2));
t116 = t110 * t112;
t114 = cos(qJ(2));
t115 = t110 * t114;
t113 = cos(qJ(3));
t111 = sin(qJ(3));
t106 = t108 * t112 * t111 - t110 * t113;
t105 = (-t107 * t116 + t109 * t114) * t111 - t113 * t118;
t104 = (t107 * t114 + t109 * t116) * t111 + t113 * t117;
t1 = [0, t118, t107 * t115 + t109 * t112, t105, t105, t105; 0, -t117, t107 * t112 - t109 * t115, t104, t104, t104; 0, t110, -t108 * t114, t106, t106, t106;];
Jg_rot  = t1;
