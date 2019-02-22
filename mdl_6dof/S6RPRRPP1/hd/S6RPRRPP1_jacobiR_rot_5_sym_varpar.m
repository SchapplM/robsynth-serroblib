% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:38
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPP1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:38:31
% EndTime: 2019-02-22 10:38:31
% DurationCPUTime: 0.04s
% Computational Cost: add. (62->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t75 = qJ(4) + pkin(10);
t71 = sin(t75);
t77 = sin(qJ(3));
t82 = t77 * t71;
t73 = cos(t75);
t81 = t77 * t73;
t78 = cos(qJ(3));
t80 = t78 * t71;
t79 = t78 * t73;
t76 = qJ(1) + pkin(9);
t74 = cos(t76);
t72 = sin(t76);
t70 = t72 * t71 + t74 * t79;
t69 = t72 * t73 - t74 * t80;
t68 = t74 * t71 - t72 * t79;
t67 = t72 * t80 + t74 * t73;
t1 = [t68, 0, -t74 * t81, t69, 0, 0; t70, 0, -t72 * t81, -t67, 0, 0; 0, 0, t79, -t82, 0, 0; t67, 0, t74 * t82, -t70, 0, 0; t69, 0, t72 * t82, t68, 0, 0; 0, 0, -t80, -t81, 0, 0; -t72 * t77, 0, t74 * t78, 0, 0, 0; t74 * t77, 0, t72 * t78, 0, 0, 0; 0, 0, t77, 0, 0, 0;];
JR_rot  = t1;
