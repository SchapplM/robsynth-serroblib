% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPPR6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:23
% EndTime: 2019-02-26 22:06:23
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->10), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t108 = sin(pkin(6));
t111 = sin(qJ(1));
t119 = t111 * t108;
t110 = sin(qJ(2));
t118 = t111 * t110;
t112 = cos(qJ(2));
t117 = t111 * t112;
t113 = cos(qJ(1));
t116 = t113 * t108;
t115 = t113 * t110;
t114 = t113 * t112;
t109 = cos(pkin(6));
t107 = qJ(3) + pkin(11);
t106 = cos(t107);
t105 = sin(t107);
t1 = [0, t119, t109 * t117 + t115, 0, 0 (-t109 * t118 + t114) * t106 + t105 * t119; 0, -t116, -t109 * t114 + t118, 0, 0 (t109 * t115 + t117) * t106 - t105 * t116; 1, t109, -t108 * t112, 0, 0, t108 * t110 * t106 + t109 * t105;];
Jg_rot  = t1;
