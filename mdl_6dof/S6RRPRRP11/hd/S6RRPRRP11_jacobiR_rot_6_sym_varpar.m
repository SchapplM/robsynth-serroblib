% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:37
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:37:13
% EndTime: 2019-02-22 11:37:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (51->13), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
t86 = sin(qJ(2));
t87 = sin(qJ(1));
t93 = t87 * t86;
t85 = qJ(4) + qJ(5);
t83 = sin(t85);
t88 = cos(qJ(2));
t82 = t88 * t83;
t84 = cos(t85);
t92 = t88 * t84;
t89 = cos(qJ(1));
t91 = t89 * t86;
t90 = t89 * t88;
t81 = -t83 * t93 + t89 * t84;
t80 = t89 * t83 + t84 * t93;
t79 = t83 * t91 + t87 * t84;
t78 = -t87 * t83 + t84 * t91;
t1 = [t81, t83 * t90, 0, t78, t78, 0; t79, t87 * t82, 0, t80, t80, 0; 0, t86 * t83, 0, -t92, -t92, 0; -t80, t84 * t90, 0, -t79, -t79, 0; t78, t87 * t92, 0, t81, t81, 0; 0, t86 * t84, 0, t82, t82, 0; -t87 * t88, -t91, 0, 0, 0, 0; t90, -t93, 0, 0, 0, 0; 0, t88, 0, 0, 0, 0;];
JR_rot  = t1;
